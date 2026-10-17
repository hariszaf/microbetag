"""
For a single, beneficiary species, extract seed complementarities with all other species
in the dataset and calculate its seed scores as node A (microbetag-seed-compls).
"""
import argparse
import json

from microbetag.cli import SEED_KEGG_MODULE


def species_scores_compls(species, conf_dic, nonseeds_dic, module_related, scores_outfile, compls_outfile):
    """
    Get scores and complements for a specific model (species).

    Note:
        Since we get all pairwise combinations, we do not care about the as_donor case for a species,
        since it is calculated when the other species is the beneficiary.
    """
    from microbetag.PhyloMint.lib.CalculateIndexes import calculate_scores, extract_complements
    from microbetag.seed_complementarity import (
        generate_fixed_pairwise_comparisons,
        get_kegg_module_related,
        kegg_module_related_intersect,
    )

    modules_ms_cpd = get_kegg_module_related(SEED_KEGG_MODULE) if module_related else None

    scores, compls = set(), {}

    species_conf = conf_dic.get(species)
    if species_conf is None:
        return None, None

    as_beneficiary, _ = generate_fixed_pairwise_comparisons(species, list(conf_dic.keys()))

    for partner in [pair[1] for pair in as_beneficiary if pair[1] != species]:

        conf, non_seed = conf_dic.get(partner), nonseeds_dic.get(partner)
        if conf is None or non_seed is None:
            continue

        SeedA, SeedB, nonSeedB = set(species_conf.keys()), set(conf.keys()), set(non_seed)

        mi_coop, mi_comp = calculate_scores(SeedA, species_conf, SeedB, nonSeedB)
        scores.add(f"{species}\t{partner}\t{mi_comp}\t{mi_coop}\n")

        B_complements_A = extract_complements(SeedA, nonSeedB)
        if module_related:
            B_complements_A = kegg_module_related_intersect(B_complements_A, modules_ms_cpd)

        # extract_complements() returns a set, which is not JSON serialisable
        compls[partner] = sorted(B_complements_A) if isinstance(B_complements_A, set) else B_complements_A

    with open(scores_outfile, "a") as f:
        f.writelines(sorted(scores))

    with open(compls_outfile, "w") as f:
        json.dump(compls, f, indent=2)

    return scores, compls


def main(argv=None):
    parser = argparse.ArgumentParser(prog="microbetag-seed-compls", description=__doc__)
    parser.add_argument("species")
    parser.add_argument("nonseeds_json")
    parser.add_argument("confidence_json")
    parser.add_argument("module_related", help="'true' to keep only KEGG module related compounds")
    parser.add_argument("scores_outfile")
    parser.add_argument("compls_outfile")
    args = parser.parse_args(argv)

    with open(args.confidence_json) as f:
        conf_dic = json.load(f)
    with open(args.nonseeds_json) as f:
        nonseeds_dic = json.load(f)

    species_scores_compls(
        args.species,
        conf_dic,
        nonseeds_dic,
        args.module_related.lower() == "true",
        args.scores_outfile,
        args.compls_outfile,
    )


if __name__ == "__main__":
    main()
