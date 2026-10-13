#!/usr/bin/env python3
"""
Annotate the proteins of a genome with RAST, for a ModelSEEDpy reconstruction (recon.py).

RAST runs on the RAST server, so this step needs an internet connection; its output can be reused
for later reconstructions.

Usage:
    rast_annotate.py --faa bin_101.faa --out bin_101-rast.tsv
    rast_annotate.py --faa bin_101.faa --check      # is RAST reachable? (annotates one protein)

Output: a tab-separated file with the protein id (as recon.py reads it) and its RAST function, one protein
per line (proteins without a function are left out).

For large genome catalogues, the requests to RAST are made robust:
- the proteins are sent in batches (--batch-size);
- a batch that fails (no connection, timeout, server error, incomplete answer) is sent again, after a
  wait that doubles each time (--retries, --wait);
- the sequences need to be proteins: RAST would "annotate" nucleotides too;
- each answer is checked: every protein sent comes back, and a genome needs at least --min-annotated of its
  proteins annotated (0.05 by default; none annotated means RAST did not work);
- the output is written only once all batches went through, so a genome that failed has no output.
The exit status is 0 on success, 1 if RAST failed (worth trying again later), 2 if the input is not right.
"""

import argparse
import os
import sys
import time

from modelseedpy import RastClient

from recon import read_proteins


class RastError(Exception):
    pass


def log(message):
    print(message, file=sys.stderr, flush=True)


def annotate_batch(client, proteins, retries, wait, timeout):
    """RAST functions of a batch of proteins ({id: function}), retrying failed requests."""
    sent = {p["id"] for p in proteins}
    for attempt in range(1, retries + 2):
        try:
            client.rpc_client.timeout = timeout
            result = client.f(proteins)
            if not result or "features" not in result[0]:
                raise RastError("RAST returned no features")
            features = result[0]["features"]
            returned = {f["id"] for f in features}
            if returned != sent:
                raise RastError(f"RAST returned {len(returned & sent)} of the {len(sent)} proteins sent")
            return {f["id"]: f["function"] for f in features if f.get("function")}
        except Exception as e:  # requests' connection errors, timeouts, HTTP and RAST server errors
            if attempt > retries:
                raise RastError(f"RAST failed after {retries + 1} attempts: {e}") from e
            delay = wait * 2 ** (attempt - 1)
            log(f"RAST request failed ({type(e).__name__}: {e}); attempt {attempt + 1} in {delay} s")
            time.sleep(delay)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--faa", required=True, help="Protein sequences of the genome (e.g. Prodigal's .faa)")
    parser.add_argument("--out", help="Output: protein id and RAST function (.tsv)")
    parser.add_argument("--check", action="store_true", help="Only check that RAST annotates the first protein")
    parser.add_argument("--batch-size", type=int, default=2000, help="Proteins per request to RAST (default: 2000)")
    parser.add_argument("--retries", type=int, default=5, help="Attempts again of a failed request (default: 5)")
    parser.add_argument("--wait", type=int, default=30, help="Seconds before the first retry; doubles each time (default: 30)")
    parser.add_argument("--timeout", type=int, default=1800, help="Seconds to wait for an answer of RAST (default: 1800)")
    parser.add_argument("--min-annotated", type=float, default=0.05,
                        help="Minimum fraction of the proteins RAST needs to annotate (default: 0.05)")
    args = parser.parse_args(argv)
    if not args.check and not args.out:
        parser.error("--out is required, unless --check")

    # The protein id is the first word of each FASTA header, made safe as in recon.py
    genome = read_proteins(args.faa)
    proteins = [{"id": f.id, "protein_translation": f.seq} for f in genome.features if f.seq]
    if not proteins:
        log(f"No protein sequences in {args.faa}")
        return 2


    client = RastClient()
    try:
        if args.check:
            annotate_batch(client, proteins[:1], args.retries, args.wait, args.timeout)
            log(f"RAST is reachable: {client.rpc_client.url}")
            return 0

        # RAST also "annotates" nucleotide sequences (A, C, G and T are amino acids too): refuse them
        nucleotides = set("ACGTUN")
        n_dna = sum(1 for p in proteins if set(p["protein_translation"].upper()) <= nucleotides)
        if n_dna > len(proteins) / 2:
            log(f"{args.faa}: {n_dna} of {len(proteins)} sequences are nucleotides; RAST needs proteins (e.g. Prodigal's .faa)")
            return 2

        functions = {}
        batches = [proteins[i:i + args.batch_size] for i in range(0, len(proteins), args.batch_size)]
        for i, batch in enumerate(batches, 1):
            functions.update(annotate_batch(client, batch, args.retries, args.wait, args.timeout))
            if len(batches) > 1:
                log(f"{args.faa}: batch {i} of {len(batches)} annotated")
    except RastError as e:
        log(f"{args.faa}: {e}")
        return 1

    fraction = len(functions) / len(proteins)
    log(f"{args.faa}: {len(functions)} of {len(proteins)} proteins annotated with RAST ({fraction:.0%})")
    if fraction < args.min_annotated:
        log(f"{args.faa}: fewer than {args.min_annotated:.0%} of the proteins annotated; is it a protein FASTA file?")
        return 2

    # Written only now, and atomically, so that a failed run leaves no partial output
    tmp = f"{args.out}.tmp"
    with open(tmp, "w") as out:
        for protein in proteins:
            if protein["id"] in functions:
                out.write(f"{protein['id']}\t{functions[protein['id']]}\n")
    os.replace(tmp, args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
