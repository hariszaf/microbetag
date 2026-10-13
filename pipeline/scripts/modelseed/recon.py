#!/usr/bin/env python3
"""
Draft genome-scale metabolic model (GEM) of a genome with ModelSEEDpy, from its proteins and their RAST
annotation (rast_annotate.py).

Usage:
    recon.py --faa bin_101.faa --rast bin_101-rast.tsv --template gram_neg --out bin_101-draft.xml

Without --rast, the proteins are annotated with RAST here (this needs an internet connection).

The template is the set of ModelSEED reactions and the biomass the model is built from:
gram_neg (default), gram_pos or core. ModelSEEDpy's genome classifier, which would choose it, cannot be
loaded with the scikit-learn versions ModelSEEDpy 0.4.2 needs, so it is a parameter, as for CarveMe.

The model is not gap-filled; gapfill.py does that. Its objective is the biomass reaction, bio1.
"""

import argparse
import os
import re
import sys

import cobra
from modelseedpy import MSBuilder, MSGenome, RastClient
from modelseedpy.core.mstemplate import MSTemplateBuilder
from modelseedpy.helpers import get_template

TEMPLATES = {"gram_neg": "template_gram_neg", "gram_pos": "template_gram_pos", "core": "template_core"}


def protein_header(header):
    """
    The protein id (first word of the FASTA header) and its description. The id is made safe for the gene rules
    of the model: characters other than letters, digits, _ and . become _ (e.g. MPB:bin_2@contig_1_1 -> MPB_bin_2_contig_1_1).
    """
    seq_id, _, desc = header.partition(" ")
    return re.sub(r"[^A-Za-z0-9_.]", "_", seq_id), desc or None


def read_proteins(faa):
    """The proteins of a genome as an MSGenome, with safe ids (protein_header)."""
    return MSGenome.from_fasta(faa, h_func=protein_header)


def load_template(name):
    """A ModelSEED template, by its short name (see TEMPLATES)."""
    if name not in TEMPLATES:
        sys.exit(f"Unknown template {name!r}: choose one of {', '.join(TEMPLATES)}")
    return MSTemplateBuilder.from_dict(get_template(TEMPLATES[name])).build()


def add_rast_annotation(genome, rast_tsv):
    """Add the RAST functions of rast_annotate.py to the genome, split into terms as RastClient does."""
    n = 0
    with open(rast_tsv) as f:
        for line in f:
            if not line.strip():
                continue
            protein_id, function = line.rstrip("\n").split("\t", 1)
            if protein_id not in genome.features:
                continue
            feature = genome.features.get_by_id(protein_id)
            for term in re.split("; | / | @", function):
                feature.add_ontology_term("RAST", term)
            n += 1
    return n


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--faa", required=True, help="Protein sequences of the genome (e.g. Prodigal's .faa)")
    parser.add_argument("--rast", help="RAST annotation of the proteins (rast_annotate.py); annotated here if missing")
    parser.add_argument("--template", default="gram_neg", choices=list(TEMPLATES), help="ModelSEED template (default: gram_neg)")
    parser.add_argument("--id", help="Model id (default: the name of the .faa file)")
    parser.add_argument("--out", required=True, help="Output: the draft model (SBML .xml)")
    args = parser.parse_args(argv)

    model_id = args.id or re.sub(r"\.(faa|fasta|fa)$", "", os.path.basename(args.faa))

    genome = read_proteins(args.faa)
    if not any(f.seq for f in genome.features):
        sys.exit(f"No protein sequences in {args.faa}")
    if args.rast:
        n = add_rast_annotation(genome, args.rast)
    else:
        RastClient().annotate_genome(genome)
        n = sum(1 for f in genome.features if f.ontology_terms.get("RAST"))
    print(f"{model_id}: {n} of {len(genome.features)} proteins with a RAST function", file=sys.stderr)

    template = load_template(args.template)
    model = MSBuilder.build_metabolic_model(
        model_id, genome,
        template=template,
        allow_all_non_grp_reactions=True,
        annotate_with_rast=False,  # done above
        gapfill_model=False,       # gapfill.py
    )
    cobra.io.write_sbml_model(model, args.out)
    print(f"{model_id}: draft model with {len(model.reactions)} reactions, {len(model.metabolites)} metabolites",
          file=sys.stderr)


if __name__ == "__main__":
    main()
