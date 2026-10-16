"""
Seed sets of the GTDB ModelSEED GEMs with the current BuildGraphNetX, one JSON per genome, resumable.

For the GTDB database: the GEMs (ModelSEED, PATRIC ids) are those of Zenodo record 23368237
(gtdb_modelseed_gems.zip); the sets, merged, are precalc/seed_compl/{confidence,nonseeds}.json.gz of
build_gtdb_precalc.py. A worker is replaced after 200 models (maxtasksperchild), as COBRApy and libSBML keep memory.

    python run_seed_sets.py <microbetag dir on Elements> <workers> [file listing the models to run, one path per line]
"""

import csv
import gzip
import json
import logging
import os
import sys
import time
import traceback
from multiprocessing import Pool
from pathlib import Path

logging.disable(logging.CRITICAL)

ROOT    = Path(sys.argv[1])
GEMS    = ROOT / "gtdb_modelseed_gems"
OUT     = ROOT / "gtdb_seed_sets"
SETS    = OUT / "sets"


def genome_id(path):
    name = path.name
    for suffix in (".PATRIC.xml.gz", ".PATRIC.xml", ".xml.gz", ".xml"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name


def run(gem):
    from microbetag.PhyloMint.lib import BuildGraphNetX as B
    gid = genome_id(gem)
    out = SETS / f"{gid}.json"
    if out.exists():
        return None
    t = time.time()
    try:
        conf, seeds, nonseeds = B.getSeedSet(B.buildDG(str(gem)), maxComponentSize=2)
        tmp = out.with_suffix(".tmp")
        tmp.write_text(json.dumps({"genome": gid, "confidence": conf, "nonseeds": sorted(nonseeds)}))
        os.replace(tmp, out)
        return ["ok", gid, f"{time.time() - t:.1f}", len(seeds), len(nonseeds), ""]
    except Exception:
        return ["failed", gid, f"{time.time() - t:.1f}", "", "", traceback.format_exc().splitlines()[-1]]


def main():
    import cobra
    cobra.Configuration().solver = "glpk"
    SETS.mkdir(parents=True, exist_ok=True)
    gems = sorted(p for p in GEMS.rglob("*") if p.is_file())
    if len(sys.argv) > 3:
        chosen = set(open(sys.argv[3]).read().split())
        gems = [g for g in gems if str(g) in chosen]
    todo = [g for g in gems if not (SETS / f"{genome_id(g)}.json").exists()]
    print(f"{len(gems)} models, {len(todo)} to do", flush=True)
    log = OUT / "log.tsv"
    new = not log.exists()
    with open(log, "a", newline="") as f, Pool(int(sys.argv[2]), maxtasksperchild=50) as pool:
        w = csv.writer(f, delimiter="\t")
        if new:
            w.writerow(["status", "genome", "seconds", "seeds", "nonseeds", "error"])
        for i, row in enumerate(pool.imap_unordered(run, todo, chunksize=4), 1):
            if row:
                w.writerow(row)
                f.flush()
            if i % 500 == 0:
                print(f"{i} / {len(todo)}", flush=True)
    print("done", flush=True)


if __name__ == "__main__":
    main()
