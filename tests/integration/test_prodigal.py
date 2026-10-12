"""ORF prediction with Prodigal."""

import os

import pytest

from microbetag.tools import run_prodigal

pytestmark = pytest.mark.requires_tools("prodigal")


def test_run_prodigal(tmp_path, data_dir):
    run_prodigal(str(data_dir / "prodigal" / "input_files" / "bin_101.fa"), "bin_101", str(tmp_path))

    faa = tmp_path / "bin_101.faa"
    assert faa.exists() and os.path.getsize(faa) > 0
    assert faa.read_text().startswith(">")
