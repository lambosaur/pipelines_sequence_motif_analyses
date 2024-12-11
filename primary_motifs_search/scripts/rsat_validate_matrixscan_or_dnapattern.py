#! /usr/bin/env python


import argparse
import sys
from pathlib import Path
from typing import Dict

import Bio.Seq
import Bio.SeqIO
import pandas as pd


def argparser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(),
    )

    parser.add_argument(
        "--input_postprocessed_table",
        help="Path to postprocessed table of matrix-scan results.",
        required=True,
        type=str,
    )

    parser.add_argument(
        "--input_fasta",
        help="Path to FASTA file of sequences scanned.",
        required=True,
        type=str,
    )

    parser.add_argument(
        "--output_flag",
        help="Path to output flag file to touch if all is validated.",
        required=True,
        type=str,
    )

    parser.add_argument(
        "--n_flanking_nts",
        help="Number of flanking nucleotides to remove from the extracted sequence (used for RSAT dna-pattern).",
        required=False,
        default=0,
        type=int,
    )

    return parser


def main():
    # Goal: extract the hits' sequences from the processed table,
    # and verify they match the segments extracted from the original fasta file,
    # using the calculated coordinates.

    parser = argparser()
    args = parser.parse_args()

    input_postprocessed_table = args.input_postprocessed_table
    input_fasta = args.input_fasta
    output_flag = args.output_flag

    ########################################################################
    # PREPARE

    # Load table and sequences
    postprocessed_table = pd.read_csv(input_postprocessed_table, sep="\t", header=0)
    sequences = {}
    for sr in Bio.SeqIO.parse(input_fasta, "fasta"):
        sequences[sr.id] = sr

    # post-processed full table should contain all info for mapping sequence
    # identifiers to genomic coordinates.
    COLUMNS_TABLE = [
        "seq_id",
        "sequence_start",
        "sequence_end",
        "sequence_strand",
    ]
    map_seqid_to_coords = (
        postprocessed_table.loc[:, COLUMNS_TABLE]
        .drop_duplicates()
        .set_index("seq_id")
        .to_dict(orient="index")
    )

    # Here: transform the previous seq-records' sequences into pandas Series
    # with nucleotides indexed with their absolute genomic coordinates.
    map_seqid_to_genomic_indexed_sequence: Dict[str, pd.Series] = {}
    for seq_id, sr in sequences.items():
        if seq_id not in map_seqid_to_coords:
            # The sequence had no hits, nothing to verify.
            continue

        seq = list(str(sr.seq))
        index = list(
            range(
                map_seqid_to_coords[seq_id]["sequence_start"],
                map_seqid_to_coords[seq_id]["sequence_end"],
            )
        )
        if map_seqid_to_coords[seq_id]["sequence_strand"] == "-":
            index = index[::-1]

        indexed_seq = pd.Series(seq, index=index)
        map_seqid_to_genomic_indexed_sequence[seq_id] = indexed_seq

    ########################################################################
    # PROCESS RESULTS

    # Now, for each hit, extract the sequence and verify it matches the expected sequence.
    for i, row in postprocessed_table.iterrows():
        expected_seq_str = row["sequence"]

        if args.n_flanking_nts > 0:
            expected_seq_str = expected_seq_str[
                args.n_flanking_nts : -args.n_flanking_nts
            ]

        extracted_seq = map_seqid_to_genomic_indexed_sequence[row["seq_id"]].loc[
            list(range(row["hit_genomic_start"], row["hit_genomic_end"]))
        ]
        if row["sequence_strand"] == "-":
            extracted_seq = extracted_seq[::-1]

        extracted_seq_str = "".join(extracted_seq.to_list())

        if row["ft_strand"] == "R":
            extracted_seq = pd.Series(
                Bio.Seq.Seq(extracted_seq_str).reverse_complement(),
                index=list(extracted_seq.index),
            )
            extracted_seq_str = "".join(extracted_seq.to_list())

        if not expected_seq_str == extracted_seq_str:
            raise ValueError(
                (
                    f"At index {i}, expected sequence {expected_seq_str} "
                    f"does not match extracted sequence {extracted_seq_str}"
                )
            )

    # Passed everything.
    Path(output_flag).touch()

    return 0


if __name__ == "__main__":
    sys.exit(main())
