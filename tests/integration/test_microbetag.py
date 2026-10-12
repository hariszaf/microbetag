"""
Complete microbetag run from the command line: 8 bins and an abundance table, with the
network inferred by FlashWeave. (config_from_net.yml needs input_files/hessler_net.csv and
input_files/bin2taxonomy.tsv, which are not in input_files.tar.gz.)
"""

import gzip
import shutil
import subprocess
import sys
import tarfile

import pytest

pytestmark = [
    pytest.mark.requires_tools("prodigal", "hmmsearch", "julia"),
    pytest.mark.requires_conda_env("mtg-phenotrex"),
]


# Seed complementarity builds GEMs with CarveMe, which needs Gurobi
def test_microbetag_run(data_copy, gurobi_license):
    work = data_copy("microbetag")
    with tarfile.open(work / "input_files.tar.gz") as tar:
        tar.extractall(work)
    ko_merged = work / "input_files" / "ko_merged.txt"
    with gzip.open(f"{ko_merged}.gz", "rb") as f_in, open(ko_merged, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)

    res = subprocess.run(
        [sys.executable, "-m", "microbetag.microbetag", "--config", str(work / "config_mtg.yml")],
        cwd=work, capture_output=True, text=True,
    )
    assert res.returncode == 0, res.stderr[-2000:]
    assert list((work / "mtg_complete_output").glob("*.cx2"))
