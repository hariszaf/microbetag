#!/usr/bin/env python3
"""
Gap-fill a draft genome-scale model (GEM) with DNNGIOR: the reactions to add are scored by its neural network.

Usage:
    gapfill.py --model bin_101-draft.xml --medium MM_glu.csv --out bin_101.xml

Drafts of gapseq and ModelSEEDpy are in the ModelSEED namespace, CarveMe's in BiGG; --db auto (default)
tells them apart by their reaction ids. The objective is the model's own (bio1, Growth, ...), unless given.

The medium is either
- a gapseq medium file (.csv: compounds, name, maxFlux), e.g. dat/media/MM_glu.csv of gapseq;
- a DNNGIOR medium file (.tsv: id, name, concentration, minflux, maxflux);
- or a 2-column tab-separated file: compound id and maximum uptake, without header;
with compound ids of the model's namespace (cpd00027 for ModelSEED, glc__D for BiGG). It is applied to the
exchange reactions of the draft before DNNGIOR reads it: the medium's compounds can be taken up, the others
not; it is also given to DNNGIOR. Without a medium, every compound with an exchange reaction can be taken up.
(DNNGIOR gap-fills on the exchange bounds of the draft; a gapseq draft has all its exchanges closed, which would
leave nothing to grow on.) The gap-filled model is saved with all its exchanges open.

Works with DNNGIOR (MGXlab/DNNGIOR, develop branch; solves with Gurobi), and with dnngior-optlang (--solver).
"""

import argparse
import csv
import inspect
import os
import re
import sys
import tempfile

import cobra

EXCHANGE_SUFFIX = {"ModelSEED": "_e0", "BiGG": "_e"}


def log(message):
    print(message, file=sys.stderr, flush=True)


def namespace(model):
    """ModelSEED if most reactions have ModelSEED ids (rxn00001_c0), otherwise BiGG."""
    ms = sum(1 for r in model.reactions if re.match(r"rxn\d{5}", r.id))
    return "ModelSEED" if ms > len(model.reactions) / 2 else "BiGG"


def read_medium(path):
    """{compound id: maximum uptake} from a gapseq, DNNGIOR or 2-column medium file."""
    with open(path, newline="") as f:
        header = f.readline()
        f.seek(0)
        if "," in header and "compounds" in header:
            return {r["compounds"].strip(): abs(float(r["maxFlux"])) for r in csv.DictReader(f)}
        if "\t" in header and "minflux" in header:
            return {r["id"].strip(): abs(float(r["minflux"])) for r in csv.DictReader(f, delimiter="\t")}
        uptakes = {}
        for line in f:
            if line.strip() and not line.startswith("#"):
                cpd, flux = line.rstrip("\n").split("\t")[:2]
                uptakes[cpd.strip()] = abs(float(flux))
        return uptakes


def check_medium(model, uptakes, db, path):
    """
    An error message if the medium does not fit the model: compound ids of the other namespace (ModelSEED's
    cpd00027 for a BiGG model, or the other way round), or most of its compounds without an exchange in the model.
    """
    ids = [re.sub(r"_e0?$", "", cpd) for cpd in uptakes]
    n_seed = sum(1 for cpd in ids if re.fullmatch(r"cpd\d{5}", cpd))
    if db == "BiGG" and n_seed > len(ids) / 2:
        return f"the medium {path} has ModelSEED compound ids (cpd...), the model BiGG ids: give a BiGG medium"
    if db == "ModelSEED" and n_seed < len(ids) / 2:
        return f"the medium {path} has no ModelSEED compound ids (cpd...), as the model has: give a ModelSEED medium"
    exchanges = {r.id for r in model.exchanges}
    found = sum(1 for cpd in ids if f"EX_{cpd}{EXCHANGE_SUFFIX[db]}" in exchanges)
    if found < len(ids) / 2:
        return (f"only {found} of the {len(ids)} compounds of the medium {path} have an exchange reaction in the model; "
                "is it the medium for these models?")
    return None


def apply_medium(model, uptakes, db):
    """
    Set the exchange bounds of the model to the medium: its compounds can be taken up (at most their uptake),
    the others not; without a medium, all can. Exchanges for medium compounds the model lacks are added.
    Returns the number of medium compounds that had no exchange in the draft.
    """
    suffix = EXCHANGE_SUFFIX[db]
    for reaction in model.exchanges:
        reaction.lower_bound = -1000.0 if uptakes is None else 0.0
        reaction.upper_bound = 1000.0
    if uptakes is None:
        return 0

    added = 0
    for cpd, uptake in uptakes.items():
        cpd = re.sub(r"_e0?$", "", cpd)
        ex_id = f"EX_{cpd}{suffix}"
        if ex_id not in model.reactions:
            met_id = f"{cpd}{suffix}"
            if met_id not in model.metabolites:
                model.add_metabolites([cobra.Metabolite(met_id, compartment=suffix.lstrip("_"))])
            reaction = cobra.Reaction(ex_id, lower_bound=0, upper_bound=1000.0)
            reaction.add_metabolites({model.metabolites.get_by_id(met_id): -1})
            model.add_reactions([reaction])
            added += 1
        model.reactions.get_by_id(ex_id).lower_bound = -uptake
    return added


def medium_for_dnngior(uptakes, db):
    """The medium as DNNGIOR's medium argument: {exchange id: bounds and metabolite}."""
    suffix = EXCHANGE_SUFFIX[db]
    medium = {}
    for cpd, uptake in uptakes.items():
        cpd = re.sub(r"_e0?$", "", cpd)
        medium[f"EX_{cpd}{suffix}"] = {
            "lower_bound": -uptake, "upper_bound": 1000.0, "metabolites": {f"{cpd}{suffix}": -1.0},
        }
    return medium


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--model", required=True, help="Draft model (SBML)")
    parser.add_argument("--medium", help="Medium file (see above); every compound can be taken up if missing")
    parser.add_argument("--db", default="auto", choices=["auto", "ModelSEED", "BiGG"], help="Namespace of the model (default: auto)")
    parser.add_argument("--objective", help="Reaction to gap-fill for (default: the model's objective)")
    parser.add_argument("--solver", help="Solver of dnngior-optlang, e.g. highs (DNNGIOR itself uses Gurobi)")
    parser.add_argument("--out", required=True, help="Output: the gap-filled model (SBML)")
    args = parser.parse_args(argv)

    import dnngior  # starts Gurobi with DNNGIOR: after the arguments are checked

    model = cobra.io.read_sbml_model(args.model)
    db = namespace(model) if args.db == "auto" else args.db
    objective = args.objective or next((r.id for r in model.reactions if r.objective_coefficient), None)
    if objective is None or objective not in model.reactions:
        log(f"{args.model}: no objective reaction; give it with --objective")
        return 2

    uptakes = read_medium(args.medium) if args.medium else None
    if uptakes is not None:
        if not uptakes:
            log(f"No compounds in the medium file {args.medium}")
            return 2
        problem = check_medium(model, uptakes, db, args.medium)
        if problem:
            log(f"{args.model}: {problem}")
            return 2
    missing = apply_medium(model, uptakes, db)
    medium_name = os.path.basename(args.medium) if args.medium else "all exchanges open"
    if missing:
        log(f"{args.model}: {missing} compounds of the medium had no exchange reaction in the draft; added")

    with tempfile.TemporaryDirectory() as tmp:
        draft = os.path.join(tmp, "draft.xml")
        cobra.io.write_sbml_model(model, draft)

        # The medium is also given to DNNGIOR itself: dnngior-optlang sets the exchange bounds from it,
        # rather than from the draft
        options = dict(objectiveName=objective, dbType=db)
        if uptakes is not None:
            options["medium"] = medium_for_dnngior(uptakes, db)
        if args.solver and "solver" in inspect.signature(dnngior.Gapfill.__init__).parameters:
            options["solver"] = args.solver
        try:
            gf = dnngior.Gapfill(draft, **options)
        except FileNotFoundError as e:
            log(f"{args.model}: this DNNGIOR has no {db} database ({e.filename}); "
                "BiGG models (CarveMe) need DNNGIOR's develop branch (containers/Dockerfile.dnngior)")
            return 2

    if not hasattr(gf, "gapfilledModel"):
        log(f"{args.model}: DNNGIOR found no gap-filling for {objective} on {medium_name}")
        return 1

    gapfilled = gf.gapfilledModel
    gapfilled.id = model.id or gapfilled.id
    gapfilled.objective = objective

    # Growth on the medium; the model is saved with all its exchanges open, as the seed sets need
    with gapfilled:
        apply_medium(gapfilled, uptakes, db)
        growth = gapfilled.slim_optimize()
    for reaction in gapfilled.exchanges:
        reaction.lower_bound, reaction.upper_bound = -1000.0, 1000.0
    cobra.io.write_sbml_model(gapfilled, args.out)
    log(f"{gapfilled.id}: {len(gf.added_reactions)} reactions added by DNNGIOR ({db}, {medium_name}); "
        f"growth on the medium {growth:.4g}")
    if not growth or growth <= 0:
        log(f"WARNING: {gapfilled.id} does not grow on {medium_name} after gap-filling")
    return 0


if __name__ == "__main__":
    sys.exit(main())
