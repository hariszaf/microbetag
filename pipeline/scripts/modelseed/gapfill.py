#!/usr/bin/env python3
"""
Gap-fill a ModelSEEDpy draft model (recon.py) on a medium, so that it grows.

Usage:
    gapfill.py --model bin_101-draft.xml --template gram_neg --medium MM_glu.csv --out bin_101.xml

The reactions are added from the same template as the draft. The medium is either
- a gapseq medium file (e.g. dat/media/MM_glu.csv of gapseq): columns compounds, name and maxFlux;
- or a tab-separated file with the ModelSEED compound id and its maximum uptake, without header.
Without --medium, the model is gap-filled on ModelSEED's complete medium, in which every compound with an
exchange reaction can be taken up.
"""

import argparse
import csv
import sys

import cobra
from modelseedpy import MSBuilder, MSMedia
from modelseedpy.fbapkg.mspackagemanager import MSPackageManager

from recon import TEMPLATES, load_template


def load_medium(path):
    """An MSMedia from a gapseq medium file (.csv) or a 2-column tab-separated file."""
    uptakes = {}
    with open(path, newline="") as f:
        first = f.readline()
        f.seek(0)
        if "," in first and "compounds" in first:
            for row in csv.DictReader(f):
                uptakes[row["compounds"].strip()] = float(row["maxFlux"])
        else:
            for line in f:
                if line.strip() and not line.startswith("#"):
                    cpd, flux = line.split("\t")[:2]
                    uptakes[cpd.strip()] = float(flux)
    if not uptakes:
        sys.exit(f"No compounds in the medium file {path}")

    # Strip a compartment suffix (cpd00027_e0 -> cpd00027), as MSMedia expects compound ids
    media = MSMedia.from_dict({cpd.split("_")[0]: flux for cpd, flux in uptakes.items()})
    media.id = media.name = path
    return media


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--model", required=True, help="Draft model (SBML) of recon.py")
    parser.add_argument("--template", default="gram_neg", choices=list(TEMPLATES), help="ModelSEED template of the draft (default: gram_neg)")
    parser.add_argument("--medium", help="Medium file (see above); the complete medium if missing")
    parser.add_argument("--target", default="bio1", help="Reaction to gap-fill for (default: the biomass, bio1)")
    parser.add_argument("--out", required=True, help="Output: the gap-filled model (SBML .xml)")
    args = parser.parse_args(argv)

    model = cobra.io.read_sbml_model(args.model)
    template = load_template(args.template)
    media = load_medium(args.medium) if args.medium else MSMedia("Complete", name="Complete")

    n_draft = len(model.reactions)
    # MSBuilder.gapfill_model adds the gap-filling reactions to the model with their direction;
    # (MSGapfill.integrate_gapfill_solution of ModelSEEDpy 0.4.2 leaves out some of them)
    model = MSBuilder.gapfill_model(model, args.target, template, media)
    model.objective = args.target

    # Growth on the medium, with only its compounds taken up; the bounds are restored after the with block,
    # so the model is saved with all its exchanges open
    with model:
        MSPackageManager.get_pkg_mgr(model).getpkg("KBaseMediaPkg").build_package(media)
        growth = model.slim_optimize()

    cobra.io.write_sbml_model(model, args.out)
    print(f"{model.id}: {len(model.reactions) - n_draft} reactions added on medium {media.name}; "
          f"growth {growth:.4g}", file=sys.stderr)
    if not growth or growth <= 0:
        print(f"WARNING: {model.id} does not grow on {media.name} after gap-filling", file=sys.stderr)


if __name__ == "__main__":
    main()
