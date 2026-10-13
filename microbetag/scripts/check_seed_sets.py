#!/usr/bin/env python3
"""
Check the seed sets microbetag infers for GEMs against the models themselves (flux balance analysis).

Usage:
    python scripts/check_seed_sets.py tests/data/seed_compl/iML1515.xml.gz tests/data/seed_compl/input_files/*.xml
    python scripts/check_seed_sets.py model.xml --json details.json

For each model, per compound (compartments merged, as the seed sets are):

- seeds, non-seeds and their overlap (none expected)
- required nutrients: on the complete medium (every uptake open), the model does not grow without them, so nothing
  in the network can replace them; each of them should be a seed (completeness)
- seeds the model can make: on the complete medium, with the seed's own uptake closed, a demand reaction for it
  carries flux; these should not be seeds (precision)
- biomass precursors reached by network expansion from the seeds (and the currency metabolites, always available)

Seeds are topological: what the network cannot make from anything, needed or not. A required nutrient the network
appears to make (a reaction produces it, but its other substrates or its stoichiometry cannot be satisfied) is a
limit of the topological definition rather than of the inference. See "How microbetag gets seed sets" in the README.
"""

import argparse
import json
import re
import sys

import cobra

from microbetag.PhyloMint.lib import BuildGraphNetX

SUFFIX = re.compile(r"_([a-z]\d?)$")


def compound(met_id):
    return SUFFIX.sub("", re.sub(r"^M_", "", met_id))


CURRENCY = {compound(c) for c in BuildGraphNetX.CURRENCY}


def complete_medium(model):
    for r in model.exchanges:
        r.lower_bound = -1000


def exchanges_by_compound(model):
    by_compound = {}
    for r in model.exchanges:
        for m in r.metabolites:
            by_compound.setdefault(compound(m.id), []).append(r)
    return by_compound


def grows(value):
    return bool(value) and value == value and value > 1e-6


def required_nutrients(model):
    """Uptakes the model cannot grow without, on the complete medium (currency metabolites left out)."""
    model = model.copy()
    complete_medium(model)
    needed = []
    for c, exchanges in sorted(exchanges_by_compound(model).items()):
        if c in CURRENCY:
            continue
        with model:
            for r in exchanges:
                r.lower_bound = 0
            if not grows(model.slim_optimize()):
                needed.append(c)
    return needed


def producible(model, seeds):
    """Seeds the model can make, with their own uptake closed (their cytosolic form, if any)."""
    model = model.copy()
    complete_medium(model)
    exchanges = exchanges_by_compound(model)
    forms = {}
    for m in model.metabolites:
        forms.setdefault(compound(m.id), []).append(m)
    made = []
    for c in sorted(seeds):
        if c in CURRENCY or c not in forms:
            continue
        met = sorted(forms[c], key=lambda m: (m.compartment not in ("c", "c0"), m.id))[0]
        with model:
            for r in exchanges.get(c, []):
                r.lower_bound = 0
            model.objective = model.add_boundary(met, type="demand")
            if grows(model.slim_optimize()):
                made.append(c)
    return made


def expansion(model, seeds):
    """Compounds reached by network expansion from the seeds, over the model's reactions in their directions."""
    skip = {r.id for r in model.boundary} | {r.id for r in model.reactions if r.objective_coefficient}
    steps = []
    for r in model.reactions:
        if r.id in skip:
            continue
        subs = {compound(m.id) for m in r.reactants} - CURRENCY
        prods = {compound(m.id) for m in r.products} - CURRENCY
        if r.upper_bound > 0:
            steps.append((subs, prods))
        if r.lower_bound < 0:
            steps.append((prods, subs))
    available = set(seeds) | CURRENCY
    changed = True
    while changed:
        changed = False
        for subs, prods in steps:
            if subs <= available and not prods <= available:
                available |= prods
                changed = True
    return available


def check(path, max_component_size=2):
    model = cobra.io.read_sbml_model(str(path))
    _, seeds, nonseeds = BuildGraphNetX.getSeedSet(
        BuildGraphNetX.buildDG(str(path)), maxComponentSize=max_component_size
    )
    seeds, nonseeds = {compound(s) for s in seeds}, {compound(n) for n in nonseeds}
    needed = required_nutrients(model)
    biomass = {compound(m.id) for r in model.reactions if r.objective_coefficient for m in r.reactants} - CURRENCY
    reached = expansion(model, seeds)
    return {
        "model": str(path),
        "seeds": sorted(seeds),
        "overlap": sorted(seeds & nonseeds),
        "required": needed,
        "required not seeds": [c for c in needed if c not in seeds],
        "seeds producible": producible(model, seeds),
        "biomass not reached": sorted(biomass - reached),
        "n non-seeds": len(nonseeds),
        "n biomass precursors": len(biomass),
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("models", nargs="+", help="GEMs (SBML, .xml or .xml.gz)")
    parser.add_argument("--max-component-size", type=int, default=2, help="as microbetag uses it (default: 2)")
    parser.add_argument("--json", help="Write the details (the compounds of each check) to this file")
    args = parser.parse_args(argv)

    cobra.Configuration().solver = "glpk"
    results = []
    print(f"{'model':40s} {'seeds':>5s} {'overlap':>7s} {'required not seeds':>18s} {'seeds producible':>16s} {'biomass reached':>15s}")
    for path in args.models:
        r = check(path, args.max_component_size)
        results.append(r)
        n_bio = r["n biomass precursors"]
        print(f"{path[-40:]:40s} {len(r['seeds']):5d} {len(r['overlap']):7d} "
              f"{len(r['required not seeds']):>8d} of {len(r['required']):<6d} "
              f"{len(r['seeds producible']):>8d} of {len(r['seeds']):<4d} "
              f"{n_bio - len(r['biomass not reached']):>8d} of {n_bio:<4d}")
        sys.stdout.flush()
    if args.json:
        with open(args.json, "w") as f:
            json.dump(results, f, indent=1)


if __name__ == "__main__":
    main()
