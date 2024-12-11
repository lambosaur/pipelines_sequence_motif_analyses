from snakemake.utils import min_version


min_version("7.6")


import os


configfile: "tests/config/filepaths.yaml"
configfile: "tests/config/rsat_matrix_scan.yaml"
configfile: "tests/config/rsat_dna_pattern.yaml"


module rsat_matrix_scan:
    snakefile:
        "smk-rules/rsat_matrix_scan.smk"
    config:
        config


use rule * from rsat_matrix_scan as rms_*


module rsat_dna_pattern:
    snakefile:
        "smk-rules/rsat_dna_pattern.smk"
    config:
        config


use rule * from rsat_dna_pattern as rdp_*


rule all:
    input:
        rules.rms_all.input,
        rules.rdp_all.input,
