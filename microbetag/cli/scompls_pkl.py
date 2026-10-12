"""Convert a seed complementarities JSON to a pickled DataFrame (microbetag-scompls-pkl)."""
import argparse
import json
import pickle

import pandas as pd


def main(argv=None):
    parser = argparse.ArgumentParser(prog="microbetag-scompls-pkl", description=__doc__)
    parser.add_argument("compls_json")
    parser.add_argument("compls_pkl")
    args = parser.parse_args(argv)

    with open(args.compls_json) as f:
        compls_dict = json.load(f)

    df = pd.DataFrame.from_dict(compls_dict)

    # Replace NaNs with [] only in the float columns
    float_cols = df.select_dtypes(include="float").columns
    df[float_cols] = df[float_cols].where(df[float_cols].notna(), [[]])

    # Attention! We need df.T, otherwise source and target are swapped
    with open(args.compls_pkl, "wb") as f:
        pickle.dump(df.T, f)


if __name__ == "__main__":
    main()
