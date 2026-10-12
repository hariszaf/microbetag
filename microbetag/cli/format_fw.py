"""Format an abundance table as FlashWeave input (microbetag-format-fw)."""
import argparse
import os
import sys

import pandas as pd


def main(argv=None):
    parser = argparse.ArgumentParser(prog="microbetag-format-fw", description=__doc__)
    parser.add_argument("abundance_file")
    parser.add_argument("-o", "--output", default="flashweave_abd_table.tsv")
    args = parser.parse_args(argv)

    # TODO (Haris Zafeiropoulos, 2025-11-10): do we need those two lines now?
    os.environ["COBRAPY_CACHE_DIR"] = "/tmp/cobrapy_cache"
    os.makedirs("/tmp/cobrapy_cache", exist_ok=True)

    # Silence import-time prints
    sys.stdout = open(os.devnull, "w")
    try:
        from microbetag.utils import detect_separator
    finally:
        sys.stdout.close()
        sys.stdout = sys.__stdout__

    print(f"Input abundance file: {args.abundance_file}", file=sys.stderr)

    delimiter        = detect_separator(args.abundance_file)
    flashweave_table = pd.read_csv(args.abundance_file, sep=delimiter).iloc[:, :-1]

    float_cols = flashweave_table.select_dtypes(include=["float64"]).columns
    for col in float_cols:
        flashweave_table[col] = flashweave_table[col].astype("int64")

    flashweave_table.iloc[:, 0] = flashweave_table.iloc[:, 0].astype(str)
    flashweave_table.to_csv(args.output, sep="\t", index=False)

    print(f"Output written to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
