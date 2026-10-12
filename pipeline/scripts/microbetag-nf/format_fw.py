#!/usr/bin/env python3
"""
Format an abundance table as FlashWeave input: drop the last column (taxonomy, or the sequences
in an amplicon table), make the counts integers and write flashweave_abd_table.tsv.

Usage: format_fw.py <abundance_file>
"""

import sys

import pandas as pd

if len(sys.argv) < 2:
    sys.exit(__doc__)

abd_file = sys.argv[1]
print(f"Input abundance file: {abd_file}", file=sys.stderr)

# sep=None lets pandas detect the delimiter (tab, comma, ...)
flashweave_table = pd.read_csv(abd_file, sep=None, engine="python").iloc[:, :-1]

float_cols = flashweave_table.select_dtypes(include=["float64"]).columns
for col in float_cols:
    flashweave_table[col] = flashweave_table[col].astype("int64")

flashweave_table.iloc[:, 0] = flashweave_table.iloc[:, 0].astype(str)
flashweave_table.to_csv("flashweave_abd_table.tsv", sep="\t", index=False)

print("Output written to flashweave_abd_table.tsv", file=sys.stderr)
