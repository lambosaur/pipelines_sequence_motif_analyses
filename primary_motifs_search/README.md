# Pipeline - Primary motif search

## Overview

This sub-directory is dedicated to scanning sequences with defined sequence
motifs and obtain their matching positions within sets of sequences.

Example of tasks:

- use motifs identified with the `primary_motifs_discovery` pipeline to scan the back the sequences used for the discovery
- scan sequences with motifs from databases
- scan with known patterns as regex or exact string match

## Content

Snakemake pipelines to run RSAT matrix-scan and dna-pattern are contained in `smk-rules`.
External scripts (called from Snakemake) or independent scripts are under `scripts`.

## Usage

The pipeline is run with Snakemake, and will generate results file for **experiments** defined in the `config/<config_file>.yaml` file, see provided template config files.

Briefly: each experiment should define a set of parameters as well as identifiers for input files, which will be used by the pipeline for generating the expected outputs.

The identifiers to input files should be keys in the `config/filepaths.yaml` file

NOTE: one can have multiple config files to separate groups of experiments for readability. The Snakefile needs to be modified accordingly by adding the `configfile: "path/to/other_config.yaml"` directive, and Snakemake will simply merge the configurations.

```bash
snakemake --cores 3 --use-conda --reason all --rerun-incomplete --printshellcmds --keep-going --latency-wait 120 --snakefile ./Snakefile ;
```

To test the pipeline:

```bash
snakemake --cores 3 --use-conda --reason all --rerun-incomplete --printshellcmds --keep-going --latency-wait 120 --snakefile ./Snakefile_test.smk ;
```

You can also run the pipeline with Singularity/Apptainer:

```bash
snakemake --cores 3 --use-singularity --reason all --rerun-incomplete --printshellcmds --keep-going --latency-wait 120 --snakefile ./Snakefile_test.smk ;
```
