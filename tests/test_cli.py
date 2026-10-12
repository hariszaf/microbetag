"""
Tests for the console scripts used by the microbetag-nf Nextflow pipeline (microbetag.cli).

Expected values come from the outputs of test_seed_compl.py, stored in test_data/test_seed_compl.
"""

import json
import os
import subprocess
import sys

import pytest

from microbetag.cli import seed_compls, seed_sets

root_dir  = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
seed_data = os.path.join(root_dir, "test_data", "test_seed_compl")

CLI_MODULES = [
    "annotate", "format_fw", "format_manta", "pcompl", "pcompl_ext",
    "seed_sets", "seed_compls", "scompls_pkl", "ssets_pkl",
]


@pytest.mark.parametrize("module", CLI_MODULES)
def test_help(module):
    """Every command imports cleanly and parses --help."""
    res = subprocess.run(
        [sys.executable, "-m", f"microbetag.cli.{module}", "--help"],
        capture_output=True, text=True,
    )
    assert res.returncode == 0, res.stderr
    assert "usage:" in res.stdout


def test_seed_sets(tmp_path):
    outfile = tmp_path / "bin_101_sets.json"
    seed_sets.main([os.path.join(seed_data, "input_files", "bin_101.xml"), "BiGG", str(outfile)])

    with open(outfile) as f:
        result = json.load(f)
    with open(os.path.join(seed_data, "output_files", "SeedSetDic.json")) as f:
        expected = json.load(f)

    assert result["base_name"] == "bin_101"
    assert set(result["seed_set"]) == set(expected["bin_101"])


def test_seed_sets_bad_namespace(tmp_path):
    with pytest.raises(SystemExit):
        seed_sets.main(["model.xml", "kegg", str(tmp_path / "out.json")])


def test_seed_compls(tmp_path):
    scores, compls = tmp_path / "scores.tsv", tmp_path / "compls.json"
    seed_compls.main([
        "bin_101",
        os.path.join(seed_data, "output_files", "nonSeedSetDic.json"),
        os.path.join(seed_data, "output_files", "confidenceDic.json"),
        "false",
        str(scores),
        str(compls),
    ])

    assert scores.read_text().splitlines() == ["bin_101\tGPB:bin_000083\t0.5\t0.22"]
    with open(compls) as f:
        assert list(json.load(f)) == ["GPB:bin_000083"]
