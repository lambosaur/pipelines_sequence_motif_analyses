#! /usr/bin/env python

import argparse
import os
import sys
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Union

import pandas as pd
import pandera as pa
from typing_extensions import Self

sys.path.insert(0, os.path.normpath(Path(__file__).parent / "../src/"))

from lib_rsat_misc import (
    Bed6Schema,
    DnaPatternProcessedSchema,
    DnaPatternRawSchema,
    HitGenomicCoordinatesSchema,
    MatrixScanProcessedSchema,
    MatrixScanRawSchema,
    load_bed6,
    parse_genomic_intervals_from_seq_id,
    read_dnapattern_raw_results,
    read_matrixscan_raw_results,
)
from pandera.typing import DataFrame


def resolve_hit_genomic_strand(rsat_strand: str, sequence_strand: str) -> str:
    """Define the final gemomic strand of a RSAT hit considering the source sequence strand.

    Define the final strand of the hit basing on the source sequence strand and the
    D or R matrix-scan result strand. e.g. a hit on the R complement strand of a sequence
    from the (-) strand of the genome will be annotated as occurring on the + strand.

    NOTE: at some point I had considered only reporting the source sequence strand,
    while maintaining the D/R hit strand in the name field of the BED6 format.
    """

    if sequence_strand not in ["+", "-"]:
        raise ValueError("Unexpected sequence strand character. (expected + or -)")
    if rsat_strand not in ["D", "R"]:
        raise ValueError("Unexpected RSAT strand character. (expected D or R)")

    if rsat_strand == "D":
        return sequence_strand
    else:
        return {"+": "-", "-": "+"}[sequence_strand]


def calculate_hit_absolute_genomic_coordinates(
    row: Mapping[str, int],
    seq_meta: Mapping[str, Union[str, int]],
) -> Dict[str, int]:
    """Expects columns "ft_start", "ft_end", "ft_strand" in row, and "start", "end" and "strand" in seq_meta.

    Return a dict with keys "genomic_start", "genomic_end" coordinates
    of a 0-based, half-open interval.
    """
    if seq_meta["strand"] == "+":
        if row["ft_strand"] == "D":
            abs_start = seq_meta["end"] + row["ft_start"]
            abs_end = seq_meta["end"] + row["ft_end"] + 1

        elif row["ft_strand"] == "R":
            abs_start = seq_meta["end"] + row["ft_start"]
            abs_end = seq_meta["end"] + row["ft_end"] + 1

        else:
            raise ValueError("Unexpected strand character in matrix-scan results.")

    elif seq_meta["strand"] == "-":
        if row["ft_strand"] == "D":
            abs_start = seq_meta["start"] - row["ft_end"] - 1
            abs_end = seq_meta["start"] - row["ft_start"]

        elif row["ft_strand"] == "R":
            abs_start = seq_meta["start"] - row["ft_end"] - 1
            # abs_start = seq_meta['start'] - (row['ft_end'] + 1)
            abs_end = seq_meta["start"] - (row["ft_start"])

        else:
            raise ValueError("Unexpected strand character in matrix-scan results.")

    else:
        raise ValueError("Unexpected strand character in sequence metadata.")

    return {"genomic_start": abs_start, "genomic_end": abs_end}


class ResultsPostprocessor(ABC):
    # This should be renamed in both dna-pattern and matrix-scan results.
    _seq_id_column: str = "seq_id"
    _bed_cols: List[str]

    @property
    def processed_table(self):
        return self._processed_table

    @property
    def seq_id_column(self):
        return self._seq_id_column

    @property
    def bed_cols(self):
        return self._bed_cols

    @property
    def raw(self):
        return self._raw

    @property
    def limits(self):
        return self._limits

    @property
    def hits(self):
        return self._hits

    @property
    def mapping_coordinates(self):
        return self._mapping_coordinates

    @classmethod
    @abstractmethod
    def from_filepath(cls, fp: os.PathLike) -> Self: ...

    @abstractmethod
    def process(self, bed: DataFrame[Bed6Schema]): ...

    def _generic_process(
        self, bed: Optional[DataFrame[Bed6Schema]] = None
    ) -> pd.DataFrame:
        """Process the raw results into a full table with genomic coordinates of hits.

        The processing is the same for both matrix-scan and dna-pattern results.

        STEPS:
        1. Get the mapping between sequence IDS and coordinates
        2. Calculate genomic coordinates of hits.
        3. Resolve the strand of the hits basing on the RSAT strand and genomic strand.
        4. Merge into a complete table.

        Args:
            bed (DataFrame[Bed6Schema]): Optional BED file to produce a mapping between
                                            unique NAME values and their genomic coordinates.
        """

        self._mapping_coordinates = self._generate_mapping_coordinates()
        hits_genomics_coordinates = (
            self._generate_table_of_absolute_coordinates_for_hits()
        )

        # Some assertions
        # TODO: ideally a way to make sure that we merge correctly the tables.
        if self.hits.shape[0] != hits_genomics_coordinates.shape[0]:
            raise ValueError("Mismatch in number of hits and genomic coordinates.")

        # Now merge everything into a final table.
        full_results_hits = pd.concat([self.hits, hits_genomics_coordinates], axis=1)

        # Merge in the parsed sequence coordinates
        indexed_parsed_seq_coordinates = pd.DataFrame(
            self.mapping_coordinates
        ).T.add_prefix("sequence_")
        # Make sure of the dtypes.
        indexed_parsed_seq_coordinates["sequence_start"] = (
            indexed_parsed_seq_coordinates["sequence_start"].astype(int)
        )
        indexed_parsed_seq_coordinates["sequence_end"] = indexed_parsed_seq_coordinates[
            "sequence_end"
        ].astype(int)

        full_results_hits = full_results_hits.merge(
            indexed_parsed_seq_coordinates,
            left_on="seq_id",
            right_index=True,
            how="left",
        )

        # Assign the hit genomic strand considering the D/R and the source sequence strand.
        full_results_hits["hit_genomic_strand"] = full_results_hits.apply(
            lambda row: resolve_hit_genomic_strand(
                row["ft_strand"], row["sequence_strand"]
            ),
            axis=1,
        )

        return full_results_hits

    def _generate_table_of_absolute_coordinates_for_hits(
        self,
    ) -> DataFrame[HitGenomicCoordinatesSchema]:
        """Calculate for each RSAT hit's relative coordinates the absolute genomic coordinates.

        This relies on the sequence ID from each row of the hits table, being
        mapped to a set of absolute genomic coordinates.

        Returns:
            DataFrame[HitGenomicCoordinatesSchema]: absolute 0-based open-ended genomic coordinates.
        """
        genomic_coordinates = pd.DataFrame(
            self.hits.apply(
                lambda row: calculate_hit_absolute_genomic_coordinates(
                    row=row,
                    seq_meta=self.mapping_coordinates[row["seq_id"]],
                ),
                axis=1,
            ).to_list()
        ).add_prefix("hit_")

        return genomic_coordinates

    def _generate_mapping_coordinates(
        self, bed: Optional[DataFrame[Bed6Schema]] = None
    ) -> Dict[str, Dict[str, int]]:
        """Generate a dictionary mapping sequence identifiers to genomic coordinates.
        If a bedfile is provided then the indexing is done on the name field,
        assuming uniquely identifiable intervals.

        If not, then the sequence identifiers are parsed from the "limits" rows
        in the results table. There, the sequence identifiers are expected to be
        formatted with genomic coordinates in the format "{chrom}:{start}-{end}:{strand}..."
        or "{chrom}:{start}-{end}({strand})...".
        """
        if bed is not None:
            # Verify uniqueness of the name, as we want to use it to map
            # the matrix-scan results to their genomic intervals.
            if bed["name"].duplicated().any():
                raise ValueError("BED file contains duplicated names.")

            # Verify that all the intervals from the results table are
            # present in the BED.
            if not set(self.limits[self.seq_id_column].unique()) <= set(
                bed["name"].unique()
            ):
                raise ValueError(
                    "BED file does not contain the same intervals as the results."
                )

            map_seq_identifier_to_genomic_coordinates = (
                bed.set_index("name")
                .loc[:, ["chrom", "start", "end", "strand"]]
                .to_dict(orient="index")
            )

        else:
            extracted_genomic_intervals = parse_genomic_intervals_from_seq_id(
                self.limits[self.seq_id_column]
            )

            map_seq_identifier_to_genomic_coordinates = (
                extracted_genomic_intervals.set_index("name")
                .loc[:, ["chrom", "start", "end", "strand"]]
                .to_dict(orient="index")
            )

        return map_seq_identifier_to_genomic_coordinates

    def export_full_table(self, output_fp: os.PathLike):
        if self.processed_table is None:
            raise ValueError("Table has not been processed yet.")

        self.processed_table.to_csv(
            output_fp,
            sep="\t",
            index=False,
            header=True,
        )

    def export_as_bed(self, output_fp: os.PathLike):
        if self.processed_table is None:
            raise ValueError("Table has not been processed yet.")

        final_cols = self.bed_cols

        # Remove potential duplicates.
        # This may happen e.g. when two motifs were identified with the same
        # ID but different patterns, which end up identifying the same hit region.
        final = (
            self.processed_table.loc[:, final_cols]
            .drop_duplicates()
            .sort_values(by=["sequence_chrom", "hit_genomic_start"])
        )

        final.to_csv(output_fp, header=False, index=False, sep="\t")


# DnaPatternResults: structure to process DNA-pattern results.
class MatrixScanResults(ResultsPostprocessor):
    _bed_cols = [
        "sequence_chrom",
        "hit_genomic_start",
        "hit_genomic_end",
        "ft_name",
        "weight",
        "hit_genomic_strand",
    ]

    def __init__(self, df: DataFrame[MatrixScanRawSchema]):
        self._raw = df

        # Extract limits, either for parsing genomic coordinates of for merging.
        self._limits = df[df["ft_type"] == "limit"].copy().reset_index(drop=True)

        # Extract motif hits
        self._hits = df[df["ft_type"] == "site"].copy().reset_index(drop=True)

        self._mapping_coordinates = None
        self._processed_table = None

    @classmethod
    def from_filepath(cls, fp: os.PathLike) -> Self:
        return cls(read_matrixscan_raw_results(fp))

    @pa.check_types(lazy=True)
    def process(
        self, bed: Optional[DataFrame[Bed6Schema]] = None
    ) -> DataFrame[MatrixScanProcessedSchema]:
        """Apply processing steps to the raw results.

        Note: this is redefined from the parent `_process` method
        so as to annotate the output type.

        Args:
            bed (Optional[DataFrame[Bed6Schema]], optional): BED dataframe to extract coordinates. Defaults to None.

        Returns:
            DataFrame[MatrixScanProcessedSchema]: full processed table ; also available at `self.processed_table`.
        """
        self._processed_table = self._generic_process(bed=bed)
        return self._processed_table


# DnaPatternResults: structure to process DNA-pattern results.
class DnaPatternResults(ResultsPostprocessor):
    _bed_cols = [
        "sequence_chrom",
        "hit_genomic_start",
        "hit_genomic_end",
        "pattern_id",
        "score",
        "hit_genomic_strand",
    ]

    def __init__(self, df: DataFrame[DnaPatternRawSchema]):
        self._raw = df

        # Extract limits, either for parsing genomic coordinates of for merging.
        self._limits = df[df["pattern_id"] == "START_END"].copy().reset_index(drop=True)

        # Extract pattern hits
        self._hits = df[df["pattern_id"] != "START_END"].copy().reset_index(drop=True)

        self._mapping_coordinates = None
        self._processed_table = None

    @classmethod
    def from_filepath(cls, fp: os.PathLike) -> Self:
        return cls(read_dnapattern_raw_results(fp))

    @pa.check_types(lazy=True)
    def process(
        self, bed: Optional[DataFrame[Bed6Schema]] = None
    ) -> DataFrame[DnaPatternProcessedSchema]:
        self._processed_table = self._generic_process(bed=bed)
        return self._processed_table


def argparser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Postprocess matrix-scan or dna-pattern results to BED format of hit coordinates.\n\n"
            #
            "The script will calculate the absolute genomic coordinates of the hits "
            "from the [-RELATIVE_START, -RELATIVE_END] in the raw table of results, "
            "using either parsed genomic coordinates of sequences from the relevant column "
            "(SeqID or #seq_id, respectively in the dna-pattern or matrix-scan results), "
            "or from the provided BED file of intervals.\n"
            #
            "When parsing genomic coordinates from the sequence IDs, the expected format "
            "of these sequence identifiers is '{chrom}:{start}-{end}:{strand}...', "
            "or '{chrom}:{start}-{end}({strand})...'.\n"
            "which are extracted by RSAT from the headers in the FASTA. "
            "**IMPORTANT**: the coordinates are expected to be **0-based, "
            "half-open intervals: [start:end), "
            "rather than 1-based, closed intervals: [start:end]**\n"
            #
            "Another assumption: negative-strand sequences are reverse-complemented. "
            "The script does not deal with non-reverse complemented sequences at the moment."
        ),
    )

    parser.add_argument(
        "--source_type",
        help="Type of source of the raw results.",
        required=True,
        type=str,
        choices=["matrixscan", "dnapattern"],
    )

    parser.add_argument(
        "--input_raw_results_fp",
        help="Path to file generated by RSAT matrix-scan or dna-pattern.",
        required=True,
        type=str,
    )

    parser.add_argument(
        "--input_bed_fp",
        help="Path to BED file of intervals for scanned sequences.",
        required=False,
        type=str,
        default=None,
    )

    parser.add_argument(
        "--output_file",
        help="Path to output file to generate.",
        required=True,
        type=str,
    )
    parser.add_argument(
        "--as_full_table",
        help="Output a complete table with parsed info, rather than BED6 format.",
        required=False,
        action="store_true",
    )

    return parser


################################################################################
# MAIN
################################################################################


def main():
    parser = argparser()
    args = parser.parse_args()

    # Load the correct structure
    if args.source_type == "matrixscan":
        structured_results = MatrixScanResults.from_filepath(args.input_raw_results_fp)

    elif args.source_type == "dnapattern":
        structured_results = DnaPatternResults.from_filepath(args.input_raw_results_fp)
    else:
        raise ValueError("Unexpected source type.")

    # Load the bedfile if provided.
    if args.input_bed_fp:
        bed = load_bed6(args.input_bed_fp)
    else:
        bed = None

    _ = structured_results.process(bed=bed)

    # Export
    if args.as_full_table:
        structured_results.export_full_table(args.output_file)
    else:
        structured_results.export_as_bed(args.output_file)

    return 0


if __name__ == "__main__":
    sys.exit(main())
