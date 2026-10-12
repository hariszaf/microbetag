import json
import pickle

import pandas as pd
from microbetag.helpers import MappingPaths
from microbetag.seed_complementarity import (
    get_kegg_module_related,
    kegg_module_related_intersect,
)

KEGG_MODULE_RELATED_SEEDS = "kegg_module_related_seeds.pkl"
KEGG_MODULE_RELATED_NONSEEDS = "kegg_module_related_nonseeds.pkl"

mappings = MappingPaths()

modelseed_compounds_of_interest = get_kegg_module_related(mappings.seed_ko_mo)

# NonSeeds =======================================


def nonseeds(all_nonseeds_file):

    with open(all_nonseeds_file, "r") as f:
        nonseeds_dic = json.load(f)

    nonseeds_dic_cl = {
        # patric: set(j.split("_")[1] for j in nonseeds)
        patric: {j.split("_")[1] if j.count("_") >= 2 else j for j in nonseeds}
        for patric, nonseeds in nonseeds_dic.items()
    }

    kegg_module_related_nonseeds = {
        patric: [
            kegg_module_related_intersect(nonseeds, modelseed_compounds_of_interest)
        ]
        for patric, nonseeds in nonseeds_dic_cl.items()
    }

    df = pd.DataFrame.from_dict(kegg_module_related_nonseeds).T
    df.index.name, df.columns = "PATRIC", ["nonSeedSet"]

    with open(KEGG_MODULE_RELATED_NONSEEDS, "wb") as f:
        pickle.dump(df, f)


# Seeds =======================================


def seeds(all_seeds_file):

    with open(all_seeds_file, "r") as f:
        seeds_dic = json.load(f)

    kegg_module_related_seeds = {}
    for patric, conf in seeds_dic.items():
        tmp_seeds = conf.keys()
        tmp_seeds = set(x.split("_")[1] if x.count("_") >= 2 else x for x in tmp_seeds)
        cl_seeds = kegg_module_related_intersect(
            tmp_seeds, modelseed_compounds_of_interest
        )
        kegg_module_related_seeds[patric] = [cl_seeds]

    df = pd.DataFrame.from_dict(kegg_module_related_seeds).T

    df.index.name, df.columns = "PATRIC", ["SeedSet"]

    with open(KEGG_MODULE_RELATED_SEEDS, "wb") as f:
        pickle.dump(df, f)
