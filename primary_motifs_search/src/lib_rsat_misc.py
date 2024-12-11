#! /usr/bin/env python

import os
import re

import pandas as pd
import pandera as pa

# import pandera.extensions as pa_extensions
from pandera.typing import DataFrame, Series

###############################################################################
# MISC


def parse_genomic_intervals_from_seq_id(seq_id_series: pd.Series) -> pd.DataFrame:
    # Sequence identifiers should be formatted so as to extract genomic coordinate of intervals.
    # Currrent expected format : "{chrom}:{start}-{end}:{strand}..."

    # "{chrom}:{start}-{end}:{strand}..."
    regex_1 = re.compile(r"(chr\w+):(\d+)-(\d+):(\+|-).*")

    # bedtools getfasta format.
    # TODO: actually it tends to put the name field BEFORE the coordinates.
    # "{chrom}:{start}-{end}({strand})..."
    regex_2 = re.compile(r"(chr\w+):(\d+)-(\d+)\((\+|-)\).*")

    # Verify uniqueness
    if seq_id_series.duplicated().any():
        raise ValueError("Sequence identifiers are not unique.")

    FOUND_REGEX = False
    for regex in [regex_1, regex_2]:
        if seq_id_series.str.match(regex).all():
            FOUND_REGEX = True
            break

    if not FOUND_REGEX:
        raise ValueError("Some sequence identifiers do not match the expected format.")

    extracted_genomic_intervals = seq_id_series.str.extract(regex)
    extracted_genomic_intervals.columns = ["chrom", "start", "end", "strand"]

    extracted_genomic_intervals["name"] = seq_id_series.values

    extracted_genomic_intervals["start"] = extracted_genomic_intervals["start"].astype(
        int
    )
    extracted_genomic_intervals["end"] = extracted_genomic_intervals["end"].astype(int)

    accepted_strand_char = ["-", "+"]
    if not ((extracted_genomic_intervals["strand"].isin(accepted_strand_char)).all()):
        raise ValueError(
            "Some parsed strand characters are not in the accepted set: [+,-]."
        )

    return extracted_genomic_intervals


class HitGenomicCoordinatesSchema:
    hit_genomic_start: Series[int] = pa.Field(ge=0)
    hit_genomic_end: Series[int] = pa.Field(gt=0)


class MiscProcessedAdditionalFieldsSchema(pa.DataFrameModel):
    sequence_chrom: Series[str]
    sequence_start: Series[int] = pa.Field(ge=0)
    sequence_end: Series[int] = pa.Field(gt=0)
    sequence_strand: Series[str] = pa.Field(isin=["+", "-"])
    hit_genomic_start: Series[int] = pa.Field(ge=0)
    hit_genomic_end: Series[int] = pa.Field(gt=0)
    hit_genomic_strand: Series[str] = pa.Field(isin=["+", "-"])


###############################################################################
# BED


class Bed6Schema(pa.DataFrameModel):
    chrom: Series[str]
    start: Series[int]
    end: Series[int]
    name: Series[str]
    score: Series
    strand: Series[str]


@pa.check_types(lazy=True)
def load_bed6(bed_fp: os.PathLike) -> DataFrame[Bed6Schema]:
    # Check if first row is header.
    bed_first_row = pd.read_csv(bed_fp, sep="\t", header=None, index_col=None, nrows=1)
    if str(bed_first_row[1]) == "start":
        header = 0
    else:
        header = None
    bed = pd.read_csv(bed_fp, sep="\t", header=header, index_col=None)
    bed.columns = ["chrom", "start", "end", "name", "score", "strand"]
    return bed


###############################################################################
# MATRIX SCAN


class MatrixScanRawSchema(pa.DataFrameModel):
    """Schema for the raw output results output from RSAT matrix-scan.

    NOTE: the names of the columns are slightly different from
    the ones in the source file:

    - '#seq_id' -> 'seq_id'
    - 'strand' -> 'ft_strand'
    - 'start' -> 'ft_start'
    - 'end' -> 'ft_end'

    You schould rename the columns accordingly when reading the file.
    """

    seq_id: Series[str]  # Renamed from '#seq_id'
    ft_type: Series[str] = pa.Field(isin=["site", "limit"])
    ft_name: Series[str]
    ft_strand: Series[str] = pa.Field(isin=["R", "D"])  # Renamed from "strand"

    # We expect start and end to be relative to the END of the sequence.
    ft_start: Series[int] = pa.Field(lt=-1)  # Renamed from "start"
    ft_end: Series[int] = pa.Field(le=-1)  # Renamed from "end"

    sequence: Series[str]
    weight: Series[float]
    Pval: Series[float]
    ln_Pval: Series[float]
    sig: Series[float] = pa.Field(nullable=True)


@pa.check_types(lazy=True)
def read_matrixscan_raw_results(
    matrixscan_fp: os.PathLike,
) -> DataFrame[MatrixScanRawSchema]:
    comment = ";"
    separator = "\t"

    df = pd.read_csv(
        matrixscan_fp, sep=separator, comment=comment, header=0, index_col=None
    )
    # Rename some columns so as to avoid conflicts with BED columns.
    MATRIX_SCAN_COLUMNS_RENAMING = {
        "#seq_id": "seq_id",
        "strand": "ft_strand",
        "start": "ft_start",
        "end": "ft_end",
    }
    df.rename(columns=MATRIX_SCAN_COLUMNS_RENAMING, inplace=True)
    return df


class MatrixScanProcessedSchema(
    MatrixScanRawSchema, MiscProcessedAdditionalFieldsSchema
):
    pass


###############################################################################
# DNA-PATTERN


class DnaPatternRawSchema(pa.DataFrameModel):
    """Schema for the raw output results output from RSAT matrix-scan.

    NOTE: the names of the columns are slightly different from
    the ones in the source file:

    - '; PatID' -> 'pattern_id'
    - 'Strand' -> 'ft_strand'
    - 'Start' -> 'ft_start'
    - 'End' -> 'ft_end'
    - 'SeqID' -> 'seq_id'
    - 'matching_seq' -> 'sequence'
    - 'Score' -> 'score'

    You schould rename the columns accordingly when reading the file.
    """

    pattern_id: Series[str]
    ft_strand: Series[str] = pa.Field(isin=["D", "R", "DR"])
    Pattern: Series[str]
    seq_id: Series[str]
    ft_start: Series[int]
    ft_end: Series[int]
    sequence: Series[str]
    score: Series[float]


class DnaPatternProcessedSchema(
    DnaPatternRawSchema, MiscProcessedAdditionalFieldsSchema
):
    pass


@pa.check_types(lazy=True)
def read_dnapattern_raw_results(
    dnapattern_fp: os.PathLike,
) -> DataFrame[DnaPatternRawSchema]:
    comment = ";"
    separator = "\t"

    # The output file contains multiple comment rows that start with ";".
    # Other rows correspond to results blocks.
    # The output file contains multiple blocks of results separated by one or
    # more comment lines. The comment line right before a block identifies
    # the type of results.
    # The block for matches should be identified by the line "; Matching positions".
    # Another block is identified by the line "; statistics  max  matched percent".
    #
    # We are only interested in the matches within the sequences.
    # NOTE: the header comes right after the line "; Matching positions".
    block_identifier_row = "; Matching positions"
    with open(dnapattern_fp, "r") as f:
        for line in f:
            if line.startswith(block_identifier_row):
                break
        # Get the next line, which is the header.
        header = next(f).strip().split(separator)
        # Read the results block.
        results_lines = []
        for line in f:
            if line.startswith(comment):
                break
            results_lines.append(line.strip().split(separator))

        df = pd.DataFrame(results_lines, columns=header)

    df.rename(
        columns={
            "; PatID": "pattern_id",
            "Start": "ft_start",
            "End": "ft_end",
            "Strand": "ft_strand",
            "SeqID": "seq_id",
            "matching_seq": "sequence",
            "Score": "score",
        },
        inplace=True,
    )

    # Apply dtypes
    df["ft_start"] = df["ft_start"].astype(int)
    df["ft_end"] = df["ft_end"].astype(int)
    df["score"] = df["score"].astype(float)

    return df
