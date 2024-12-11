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


PIPELINE_NAME = "rsat_dna_pattern"


include: "utils.smk"


###############################################################################
# CONFIGURATION AND PARSING OF CONFIG FILE


@dataclass
class DnaPatternParamsDefaults:
    """This class holds the default values for the parameters of the dna-pattern tool.

    These parameters are parsed from the config, expected under the key chain
    `rsat_dna_pattern.defaults`.

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
    mask: str = glom(
        target=config,
        spec=f"{PIPELINE_NAME}.defaults.mask",
        default="non-dna",
    )

    merge_overlapping: bool = glom(
        target=config,
        spec=f"{PIPELINE_NAME}.defaults.merge_overlapping",
        default=False,
    )

    return_n_flanking_nt: int = glom(
        target=config,
        spec=f"{PIPELINE_NAME}.defaults.return_n_flanking_nt",
        default=0,
    )


DNA_PATTERN_PARAMS_DEFAULTS = DnaPatternParamsDefaults()


class DnaPatternParams(BaseModel):
    both_strands: Optional[bool] = DNA_PATTERN_PARAMS_DEFAULTS.both_strands
    mask: Optional[str] = DNA_PATTERN_PARAMS_DEFAULTS.mask
    merge_overlapping: Optional[bool] = DNA_PATTERN_PARAMS_DEFAULTS.merge_overlapping
    return_n_flanking_nt: Optional[int] = DNA_PATTERN_PARAMS_DEFAULTS.return_n_flanking_nt


class DnaPatternExperiment(BaseModel):
    description: str

    sequences_filepath_id: str
    pattern_id: str
    # patterns here is either a regex (which will be labeled with the pattern_id in the results)
    # or empty / not provided. If not provided, then `pattern_id` should be
    # a key in the `filepaths` config that points to a file containing the patterns.
    pattern: Optional[str] = None

    params: DnaPatternParams


class ConfigDnaPatternExperiments(BaseModel):
    # Each run is identified by a unique key in the config file.
    experiments: Dict[str, DnaPatternExperiment]


# Define a model for the overall configuration of experiments to run.
CONFIG_EXPERIMENTS = ConfigDnaPatternExperiments(
    **{
        "experiments": config[PIPELINE_NAME]["experiments"]
        if PIPELINE_NAME in config
        else {}
    }
)


###############################################################################
# DEFINITIONS


def format_dnapattern_commandline_pattern_arguments(wildcards: Wildcards) -> str:
    """Format either a "-pl <pattern_filepath>" or "-p <pattern_string> -id <pattern_id>" depending on config."""
    pattern_id = CONFIG_EXPERIMENTS.experiments[wildcards.EXPERIMENT_ID].pattern_id
    pattern = CONFIG_EXPERIMENTS.experiments[wildcards.EXPERIMENT_ID].pattern

    if pattern:
        return f"-p {pattern} -id {pattern_id}"
    else:
        pattern_fp = prefix_with_basedir(config["filepaths"][pattern_id])
        return f"-pl {pattern_fp}"


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


##
rule all:
    input:
        bed_hits=expand(
            make_output_filepath("{EXPERIMENT_ID}/hits.bed"),
            EXPERIMENT_ID=config["rsat_dna_pattern"]["experiments"].keys(),
        ),
        validation_flags=expand(
            make_output_filepath("{EXPERIMENT_ID}/dna_pattern_results.validated.flag"),
            EXPERIMENT_ID=config["rsat_dna_pattern"]["experiments"].keys(),
        ),


##
rule dna_pattern:
    singularity:
        RSAT_DOCKERHUB_IMAGE
    conda:
        RSAT_CONDA_ENV
    log:
        LOGSDIR_EXEC / "{EXPERIMENT_ID}.rsat_dna_pattern.log",
    threads: 1
    input:
        fasta=lambda wildcards: get_mandatory_filepath_from_experiments_config(
            wildcards=wildcards,
            keychain=f"{wildcards.EXPERIMENT_ID}.sequences_filepath_id",
        ),
        optional_pattern_file=lambda wildcards: get_optional_filepath_from_experiments_config(
            wildcards=wildcards, keychain=f"{wildcards.EXPERIMENT_ID}.pattern_id"
        ),
    output:
        fp=make_output_filepath("{EXPERIMENT_ID}/dna_pattern_results.raw.txt"),
    params:
        pattern_args=lambda wildcards: format_dnapattern_commandline_pattern_arguments(
            wildcards
        ),
        #
        both_strands=lambda wildcards: ["-1str", "-2str"][
            int(
                CONFIG_EXPERIMENTS.experiments[
                    wildcards.EXPERIMENT_ID
                ].params.both_strands
            )
        ],
        #
        mask=lambda wildcards: CONFIG_EXPERIMENTS.experiments[
            wildcards.EXPERIMENT_ID
        ].params.mask,
        #
        return_fields="sites,limits,stats",
        #
        merge_overlapping=lambda wildcards: CONFIG_EXPERIMENTS.experiments[
            wildcards.EXPERIMENT_ID
        ].params.merge_overlapping,
        #
        return_n_flanking_nt=lambda wildcards: CONFIG_EXPERIMENTS.experiments[
            wildcards.EXPERIMENT_ID
        ].params.return_n_flanking_nt,
    shell:
        """
        rsat dna-pattern \
            -i {input.fasta} \
            -format fasta \
            \
            {params.pattern_args} \
            \
            -o {output.fp} \
            \
            {params.both_strands} \
            -return {params.return_fields} \
            -match_format table \
            -N {params.return_n_flanking_nt} \
            -origin "-0" \
            \
            -v \
            1>{log} 2>&1 \
            ;
        """


#
## rule postprocess_matrixscan_results: parse the matrix-scan table into a BED6 file of motif hits (or full-table if required)
rule postprocess_dnapattern_results:
    threads: 1
    input:
        dna_pattern_results_fp=lambda wildcards: make_output_filepath(
            f"{wildcards.EXPERIMENT_ID}/dna_pattern_results.raw.txt"
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
        table=make_output_filepath("{EXPERIMENT_ID}/dna_pattern_results.parsed.txt"),
    params:
        option_bed_fp=lambda wildcards, input: (
            "--input_bed_fp " + input.optional_bed
            if isinstance(input.optional_bed, str)
            else ""
        ),
    shell:
        """
        python scripts/rsat_postprocess_matrixscan_or_dnapattern.py \
            --source_type "dnapattern" \
            --input_raw_results_fp {input.dna_pattern_results_fp} \
            {params.option_bed_fp} \
            --output_file {output.bed} \
            ;
        python scripts/rsat_postprocess_matrixscan_or_dnapattern.py \
            --source_type "dnapattern" \
            --input_raw_results_fp {input.dna_pattern_results_fp} \
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
            "{EXPERIMENT_ID}/dna_pattern_results.parsed.txt"
        ),
        fasta=lambda wildcards: get_mandatory_filepath_from_experiments_config(
            wildcards=wildcards,
            keychain=f"{wildcards.EXPERIMENT_ID}.sequences_filepath_id",
        ),
    output:
        flag=make_output_filepath("{EXPERIMENT_ID}/dna_pattern_results.validated.flag"),
    params:
        return_n_flanking_nt=lambda wildcards: CONFIG_EXPERIMENTS.experiments[
            wildcards.EXPERIMENT_ID
        ].params.return_n_flanking_nt,
    shell:
        """
        python scripts/rsat_validate_matrixscan_or_dnapattern.py \
            --input_postprocessed_table {input.postprocessed_table} \
            --input_fasta {input.fasta} \
            --output_flag {output.flag} \
            --n_flanking_nts {params.return_n_flanking_nt} \
            ;
        """


#
