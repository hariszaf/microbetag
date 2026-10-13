#!/usr/bin/env python3
"""
Annotate the proteins of a genome with RAST, for a ModelSEEDpy reconstruction (recon.py).

RAST runs on the RAST server, so this step needs an internet connection; its output can be reused
for later reconstructions.

Usage:
    rast_annotate.py --faa bin_101.faa --out bin_101-rast.tsv

Output: a tab-separated file with the protein id (as recon.py reads it) and its RAST function, one protein
per line (proteins without a function are left out).
"""

import argparse
import sys

from modelseedpy import RastClient

from recon import read_proteins


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--faa", required=True, help="Protein sequences of the genome (e.g. Prodigal's .faa)")
    parser.add_argument("--out", required=True, help="Output: protein id and RAST function (.tsv)")
    args = parser.parse_args(argv)

    # The protein id is the first word of each FASTA header, made safe as in recon.py
    genome = read_proteins(args.faa)
    proteins = [{"id": f.id, "protein_translation": f.seq} for f in genome.features if f.seq]
    if not proteins:
        sys.exit(f"No protein sequences in {args.faa}")

    # Without splitting the functions into terms, so that recon.py gets RAST's whole function
    result = RastClient().f(proteins)

    n = 0
    with open(args.out, "w") as out:
        for feature in result[0]["features"]:
            function = feature.get("function")
            if function:
                out.write(f"{feature['id']}\t{function}\n")
                n += 1

    print(f"{args.faa}: {n} of {len(proteins)} proteins annotated with RAST", file=sys.stderr)


if __name__ == "__main__":
    main()
