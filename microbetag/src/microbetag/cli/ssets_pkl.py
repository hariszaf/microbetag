"""
Intersect seed / non-seed sets with KEGG module related compounds and pickle them
(microbetag-ssets-pkl).
"""
import argparse
import json
import pickle

import pandas as pd

from microbetag.cli import SEED_KEGG_MODULE


def main(argv=None):
    parser = argparse.ArgumentParser(prog="microbetag-ssets-pkl", description=__doc__)
    parser.add_argument("--seeds", default="seeds.json")
    parser.add_argument("--nonseeds", default="nonseeds.json")
    parser.add_argument("--seeds-pkl", default="keggm_seeds.pkl")
    parser.add_argument("--nonseeds-pkl", default="keggm_nonseeds.pkl")
    args = parser.parse_args(argv)

    from microbetag.seed_complementarity import get_kegg_module_related, kegg_module_related_intersect

    modules_ms_cpd = get_kegg_module_related(SEED_KEGG_MODULE)

    for json_file, pickle_file in [(args.nonseeds, args.nonseeds_pkl), (args.seeds, args.seeds_pkl)]:
        with open(json_file) as f:
            sets = json.load(f)

        df = pd.DataFrame.from_dict(
            {k: [kegg_module_related_intersect(v, modules_ms_cpd)] for k, v in sets.items()}
        )
        with open(pickle_file, "wb") as f:
            pickle.dump(df.T, f)


if __name__ == "__main__":
    main()
