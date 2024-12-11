from snakemake.utils import min_version


min_version("7.6")


from dataclasses import dataclass
from snakemake.io import Wildcards
from pathlib import Path
from pydantic import BaseModel
from glom import glom
from typing import Dict, Optional, Union, List
from typing_extensions import Never
import datetime

PIPELINE_NAME = "rsat_matrix_scan"


include: "utils.smk"


###############################################################################
# CONFIGURATION AND PARSING OF CONFIG FILE


@dataclass
class MatrixScanParamsDefaults:
    """This class holds the default values for the parameters of the matrix-scan tool.

    These parameters are parsed from the config, expected under the key chain
    `rsat_matrix_scan.defaults`.

    NOTE: if the key chain is not found in the config, there are also defaults
    defined in this class.

    These defaults are further used by the Pydantic model `MatrixScanParams`
    that will parse the config dict for each experiment and use defaults for
    missing values.
    """

    both_strands: bool = glom(
        target=config,
        spec=f"{PIPELINE_NAME}.defaults.both_strands",
        default=True,
    )
    lower_threshold_weight: float = glom(
        target=config,
        spec=f"{PIPELINE_NAME}.defaults.lower_threshold_weight",
        default=0,
    )
    upper_threshold_pval: float = glom(
        target=config,
        spec=f"{PIPELINE_NAME}.defaults.upper_threshold_pval",
        default=0.1,
    )


MATRIX_SCAN_PARAMS_DEFAULTS = MatrixScanParamsDefaults()


class MatrixScanParams(BaseModel):
    both_strands: Optional[bool] = MATRIX_SCAN_PARAMS_DEFAULTS.both_strands
    lower_threshold_weight: Optional[
        float
    ] = MATRIX_SCAN_PARAMS_DEFAULTS.lower_threshold_weight
    upper_threshold_pval: Optional[float] = MATRIX_SCAN_PARAMS_DEFAULTS.upper_threshold_pval


class MatrixScanExperiment(BaseModel):
    description: str

    sequences_filepath_id: str
    transfac_motifs_id: str
    bed_filepath_id: Optional[str] = None

    params: MatrixScanParams


class ConfigMatrixScanExperiments(BaseModel):
    experiments: Dict[str, MatrixScanExperiment]  # Each run is identified by a unique key


# Define a model for the overall configuration of experiments to run.
CONFIG_EXPERIMENTS = ConfigMatrixScanExperiments(
    **{
        "experiments": config[PIPELINE_NAME]["experiments"]
        if PIPELINE_NAME in config
        else {}
    }
)


###############################################################################
# MISC

LOGSDIR = Path(
    prefix_with_basedir(glom(target=config, spec="filepaths.logs_dir", default="logs"))
)
LOGSDIR_EXEC = LOGSDIR / datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

RSAT_DOCKERHUB_IMAGE = "docker://biocontainers/rsat:20240828_cv1"
RSAT_SINGULARITY_ENV = (
    prefix_with_basedir(config["filepaths"]["singularity_env_rsat"])
    if "singularity_env_rsat" in config["filepaths"]
    else RSAT_DOCKERHUB_IMAGE
)
RSAT_CONDA_ENV = prefix_with_basedir(config["filepaths"]["conda_env_rsat"])


###############################################################################
# RULES


rule all:
    input:
        bed_hits=expand(
            make_output_filepath("{EXPERIMENT_ID}/hits.bed"),
            EXPERIMENT_ID=config[PIPELINE_NAME]["experiments"].keys(),
        ),
        validation_flags=expand(
            make_output_filepath("{EXPERIMENT_ID}/matrix_scan_results.validated.flag"),
            EXPERIMENT_ID=config[PIPELINE_NAME]["experiments"].keys(),
        ),


# TODO: enable different options for the background?
rule matrix_scan:
    threads: 1
    singularity:
        RSAT_DOCKERHUB_IMAGE
    conda:
        RSAT_CONDA_ENV
    log:
        LOGSDIR_EXEC / f"{PIPELINE_NAME}.{{EXPERIMENT_ID}}.matrix_scan.log",
    input:
        fasta=lambda wildcards: get_mandatory_filepath_from_experiments_config(
            wildcards=wildcards,
            keychain=f"{wildcards.EXPERIMENT_ID}.sequences_filepath_id",
        ),
        matrix_fp=lambda wildcards: get_mandatory_filepath_from_experiments_config(
            wildcards=wildcards,
            keychain=f"{wildcards.EXPERIMENT_ID}.transfac_motifs_id",
        ),
    output:
        fp=make_output_filepath("{EXPERIMENT_ID}/matrix_scan_results.raw.txt"),
    params:
        both_strands=lambda wildcards: ["-1str", "-2str"][
            int(
                CONFIG_EXPERIMENTS.experiments[
                    wildcards.EXPERIMENT_ID
                ].params.both_strands
            )
        ],
        pseudocount=1,
        decimals=1,
        upper_threshold_pval=(
            lambda wildcards: CONFIG_EXPERIMENTS.experiments[
                wildcards.EXPERIMENT_ID
            ].params.upper_threshold_pval
        ),
        lower_threshold_weight=(
            lambda wildcards: CONFIG_EXPERIMENTS.experiments[
                wildcards.EXPERIMENT_ID
            ].params.lower_threshold_weight
        ),
    shell:
        """
        rsat matrix-scan \
            -v 1 \
            -quick \
            -i {input.fasta} \
            -seq_format fasta \
            -m {input.matrix_fp} \
            -matrix_format transfac \
            -ac_as_name \
            -o {output.fp} \
            {params.both_strands} \
            -pseudo {params.pseudocount} \
            -decimals {params.decimals} \
            -n "score" \
            \
            -bginput \
            -markov 1 \
            -bg_pseudo 0.01 \
            \
            -origin "end" \
            -return limits \
            -return sites \
            \
            -return pval \
            -uth pval {params.upper_threshold_pval} \
            -lth score {params.lower_threshold_weight} \
            1> {log} 2>&1 \
            ;
        """


#
## rule postprocess_matrixscan_results: parse the matrix-scan table into a BED6 file of motif hits (or full-table if required)
rule postprocess_matrixscan_results:
    threads: 1
    input:
        matrix_scan_results_fp=lambda wildcards: make_output_filepath(
            f"{wildcards.EXPERIMENT_ID}/matrix_scan_results.raw.txt"
        ),
        fasta=lambda wildcards: get_mandatory_filepath_from_experiments_config(
            wildcards=wildcards,
            keychain=f"{wildcards.EXPERIMENT_ID}.sequences_filepath_id",
        ),
        optional_bed=lambda wildcards: get_optional_filepath_from_experiments_config(
            wildcards=wildcards, keychain=f"{wildcards.EXPERIMENT_ID}.bed_filepath_id"
        ),
    output:
        bed=make_output_filepath("{EXPERIMENT_ID}/hits.bed"),
        table=make_output_filepath("{EXPERIMENT_ID}/matrix_scan_results.parsed.txt"),
    params:
        option_bed_fp=lambda wildcards, input: (
            "--input_bed_fp " + input.optional_bed
            if isinstance(input.optional_bed, str)
            else ""
        ),
    shell:
        """
        python scripts/rsat_postprocess_matrixscan_or_dnapattern.py \
            --source_type "matrixscan" \
            --input_raw_results_fp {input.matrix_scan_results_fp} \
            {params.option_bed_fp} \
            --output_file {output.bed} \
            ;
        python scripts/rsat_postprocess_matrixscan_or_dnapattern.py \
            --source_type "matrixscan" \
            --input_raw_results_fp {input.matrix_scan_results_fp} \
            {params.option_bed_fp} \
            --output_file {output.table} \
            --as_full_table \
            ;
        """


#
## rule validate_postprocessed_results: using the full table of parsed matrix-scan results,
## check the correctness of extracted hit coordinates through the expected vs observed sequence.
rule validate_postprocessed_results:
    threads: 1
    input:
        postprocessed_table=make_output_filepath(
            "{EXPERIMENT_ID}/matrix_scan_results.parsed.txt"
        ),
        fasta=lambda wildcards: get_mandatory_filepath_from_experiments_config(
            wildcards=wildcards,
            keychain=f"{wildcards.EXPERIMENT_ID}.sequences_filepath_id",
        ),
    output:
        flag=make_output_filepath("{EXPERIMENT_ID}/matrix_scan_results.validated.flag"),
    shell:
        """
        python scripts/rsat_validate_matrixscan_or_dnapattern.py \
            --input_postprocessed_table {input.postprocessed_table} \
            --input_fasta {input.fasta} \
            --output_flag {output.flag} \
            ;
        """


#
