#!/usr/bin/env python3
"""
assembly_stats_libs.py

Compute assembly stats (num_contigs, total_length, N50, total_length_binned)
using Biopython, pandas, and Click for a clean CLI.
"""

import click
import pandas as pd
from Bio import SeqIO


def compute_n50(lengths):
    """Return the N50 of a list of lengths."""
    sorted_lens = sorted(lengths, reverse=True)
    half = sum(sorted_lens) / 2
    cum = 0
    for L in sorted_lens:
        cum += L
        if cum >= half:
            return L
    return 0


@click.command()
@click.option(
    "-i",
    "--input-fasta",
    "fasta_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="Input contigs FASTA",
)
@click.option(
    "-b",
    "--bin-file",
    "bin_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="Two-column file (contig<tab>bin) without header",
)
@click.option(
    "-o",
    "--output-csv",
    "out_csv",
    required=True,
    type=click.Path(writable=True, dir_okay=False),
    help="Output CSV (single row)",
)
def main(fasta_path, bin_path, out_csv):
    # 1) Parse contig lengths via Biopython
    records = list(SeqIO.parse(fasta_path, "fasta"))
    lengths = [len(rec.seq) for rec in records]
    contig_names = [rec.id for rec in records]

    # Build a DataFrame of contig → length
    df = pd.DataFrame({"contig": contig_names, "length": lengths})

    # 2) Read bin assignments
    df_bins = pd.read_csv(
        bin_path,
        sep=None,  # auto-detect whitespace
        engine="python",
        names=["contig", "bin"],
        usecols=["contig"],
    )

    # Merge and compute total binned length
    df_merged = pd.merge(df, df_bins, on="contig", how="inner")
    total_length_binned = df_merged["length"].sum()

    # 3) Compute summary stats
    stats = {
        "num_contigs": len(df),
        "total_length": df["length"].sum(),
        "N50": compute_n50(df["length"].tolist()),
        "total_length_binned": total_length_binned,
    }

    # 4) Write a single‐row CSV
    out_df = pd.DataFrame([stats])
    out_df.to_csv(out_csv, index=False)

    click.echo(f"✔ Stats written to {out_csv}")


if __name__ == "__main__":
    main()

