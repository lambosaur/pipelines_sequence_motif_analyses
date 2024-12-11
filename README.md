# Pipelines for sequence motif analyses

Unofficial implementation of a set of Snakemake pipelines for motifs-based analyses, notably using RSAT.

## Overview

This repository contains a set of Snakemake pipelines for sequence motif analyses, notably using RSAT:

- [Motif search in sequences](./primary_motifs_search/README.md)


## How to use

The current suggested approach is to fork the repository for each project where you want to use the pipelines.

The forked repository can then be used as a module in your project, along other analyses, and integrate project-specific changes.

## Important note on genomic coordinates

Across the pipelines, we expect sequence intervals to be defined in 0-based, open-end coordinates, and identify unique intervals through this coordinate system:

- In FASTA files, the headers should be formatted with `>{chrom}:{start}-{end}:{strand}` or `>{chrom}:{start}-{end}({strand})` of **0-based, open-end coordinates**.
- This is the default format produced by `bedtools getfasta` **without** the `-name` option, taking as input a BED file that contains 0-based coordinates.
- **CAUTION**: some sequence-manipulation libraries expect coordinates in 1-based, closed-end intervals. Be warry of this when manually generating sequences and manipulating intervals.

## Set-up

TODO: describe precisely + have both the conda and the singularity approaches.

RSAT-env:

<https://bioconda.github.io/recipes/rsat-core/README.html>

<https://hub.docker.com/r/biocontainers/rsat/tags>

<https://github.com/rsa-tools/rsat-code/tree/master/docker>

<https://html-preview.github.io/?url=https://github.com/rsa-tools/installing-RSAT/blob/master/RSAT-Docker/RSAT-Docker-tuto.html>

<https://rsa-tools.github.io/installing-RSAT/unix-install-rsat/installing_RSAT_procedure.html>

<https://github.com/rsa-tools/installing-RSAT/blob/master/RSAT-Docker/RSAT-Docker-tuto.html>

