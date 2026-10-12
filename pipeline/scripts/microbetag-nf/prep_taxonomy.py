#!/usr/bin/env python3
"""
Helpers for the taxonomy assignment of an amplicon abundance table (sequence ids in the first
column, samples, the sequences in the last column).

  prep_taxonomy.py fasta <abundance_table> <out.fasta>
      Write the sequences of the table as FASTA.

  prep_taxonomy.py merge <abundance_table> <assignments.tsv> <out.tsv>
      Replace the sequence column with the taxonomy from IDTAXA, as microbetag expects it
      (taxonomy in the last column).
"""

import sys

import pandas as pd


def read_table(path):
    # sep=None lets pandas detect the delimiter (tab, comma, ...)
    return pd.read_csv(path, sep=None, engine="python")


def to_fasta(table, out_fasta):
    df = read_table(table)
    seq_id, seq = df.columns[0], df.columns[-1]
    with open(out_fasta, "w") as f:
        for sid, s in zip(df[seq_id].astype(str), df[seq]):
            f.write(f">{sid}\n{s}\n")


def merge(table, assignments, out_tsv):
    df = read_table(table)
    seq_id, seq = df.columns[0], df.columns[-1]
    tax = pd.read_csv(assignments, sep="\t", dtype=str)
    df[seq_id] = df[seq_id].astype(str)
    merged = df.drop(columns=seq).merge(tax, left_on=seq_id, right_on="seqid", how="left")
    merged = merged.drop(columns="seqid")
    merged.to_csv(out_tsv, sep="\t", index=False)


if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1] not in ("fasta", "merge"):
        sys.exit(__doc__)
    if sys.argv[1] == "fasta" and len(sys.argv) == 4:
        to_fasta(sys.argv[2], sys.argv[3])
    elif sys.argv[1] == "merge" and len(sys.argv) == 5:
        merge(sys.argv[2], sys.argv[3], sys.argv[4])
    else:
        sys.exit(__doc__)
