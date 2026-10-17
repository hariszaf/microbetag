"""
Compute the seed set, non-seed set and seed-set confidence of an SBML model
and write them to a JSON file (microbetag-seed-sets).
"""
import argparse
import json
import sys

from microbetag.cli import METANETX_COMPOUNDS


def main(argv=None):
    parser = argparse.ArgumentParser(prog="microbetag-seed-sets", description=__doc__)
    parser.add_argument("sbml_path")
    parser.add_argument("namespace", choices=["BiGG", "modelseed"])
    parser.add_argument("outfile")
    args = parser.parse_args(argv)

    from microbetag.seed_complementarity import Ixes, bigg_to_seed_mapping_df, process_sbml

    is_bigg = args.namespace == "BiGG"
    ixes    = Ixes("M", "e" if is_bigg else "e0", "c" if is_bigg else "c0")

    bigg2seed = bigg_to_seed_mapping_df(METANETX_COMPOUNDS)
    result    = process_sbml((args.sbml_path, args.namespace, is_bigg, bigg2seed, ixes))
    if result is None:
        sys.exit(f"Failed to compute seed sets for {args.sbml_path}")

    base_name, seed_set, non_seed_set, seed_set_confidence = result
    with open(args.outfile, "w") as f:
        json.dump(
            {
                "base_name"          : base_name,
                "seed_set"           : seed_set,
                "non_seed_set"       : non_seed_set,
                "seed_set_confidence": seed_set_confidence,
            },
            f,
            indent=2,
        )


if __name__ == "__main__":
    main()
