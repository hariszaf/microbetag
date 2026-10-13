"""
Complete microbetag run from the command line on the Hessler data set: 7 bins, their network and
the outputs of the external tools (KO annotation, GEMs, phenotrex, FAPROTAX, manta), as
microbetag-nf produces them. microbetag computes the pathway and seed complementarities and
builds the annotated network.
"""

import json
import subprocess
import sys

BINS = {"bin_101", "bin_151", "bin_19", "bin_38", "bin_41", "bin_45", "bin_48"}


def test_microbetag_run(data_copy):
    work = data_copy("microbetag")

    res = subprocess.run(
        [sys.executable, "-m", "microbetag.microbetag", "--config", str(work / "config.yml")],
        cwd=work, capture_output=True, text=True,
    )
    assert res.returncode == 0, res.stderr[-3000:]

    out = work / "output"

    # Pathway complementarities for every bin
    with open(out / "pathway_complementarity" / "pathCompls.json") as f:
        assert set(json.load(f)) == BINS

    # Seed sets of every GEM, and scores for pairs of them
    with open(out / "seeds_complementarity" / "SeedSetDic.json") as f:
        assert set(json.load(f)) == BINS
    assert (out / "seeds_complementarity" / "seed_scores.tsv").stat().st_size > 0

    # The annotated network
    cx2 = list(out.glob("*.cx2"))
    assert len(cx2) == 1
    with open(cx2[0]) as f:
        net = json.load(f)
    aspects = {key: value for aspect in net for key, value in aspect.items()}
    nodes, edges = aspects["nodes"], aspects["edges"]
    assert {n["v"]["name"] for n in nodes} >= BINS

    # Each annotation layer is present: FAPROTAX, phenotrex and manta on the nodes ...
    node_attrs = {key for n in nodes if n["v"]["name"] in BINS for key in n["v"]}
    for prefix in ("faprotax::", "phendb::", "manta::cluster"):
        assert any(key.startswith(prefix) for key in node_attrs), prefix

    # ... and pathway and seed complementarities on the edges
    edge_attrs = {key for e in edges for key in e["v"]}
    for prefix in ("compl::", "seedCompl::", "seed::cooperation", "seed::competition"):
        assert any(key.startswith(prefix) for key in edge_attrs), prefix
