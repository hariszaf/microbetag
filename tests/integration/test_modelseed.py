"""
GEM reconstruction with ModelSEEDpy and gap-filling with DNNGIOR, each in its own conda environment:

- from RAST-annotated proteins (.faa)
- from genomes (.fa), which also needs RASTtk for the annotation
"""

import os
import tarfile
from pathlib import Path

import pytest

from microbetag.genres import GEMSReconstruction

pytestmark = [
    pytest.mark.requires_conda_env("mtg-modelseed"),
    pytest.mark.requires_conda_env("mtg-dnngior"),
]


def extract(archive: Path, dest: Path) -> Path:
    with tarfile.open(archive) as tar:
        tar.extractall(dest)
    return dest


class Config:
    threads       = 2
    gapfill_model = True
    gapfill_media = None


def test_modelseed_from_faa(tmp_path, data_dir, gurobi_license):
    # RAST annotations (.faa, .gto) of bin_6
    extract(data_dir / "modelseed" / "output_files.tar.gz", tmp_path)

    config = Config()
    config.sc_input_type       = "proteins_faa"
    config.for_reconstructions = str(tmp_path / "output_files" / "reconstructions")
    config.genres              = str(tmp_path / "GENREs")
    os.makedirs(config.genres)

    GEMSReconstruction(config).modelseed_reconstructions()

    assert (tmp_path / "GENREs" / "bin_6.xml").exists()


@pytest.mark.requires_tools("rast-create-genome")
def test_modelseed_from_fasta(tmp_path, data_dir, gurobi_license):
    extract(data_dir / "modelseed" / "input_files.tar.gz", tmp_path)

    config = Config()
    config.sc_input_type   = "bins_fasta"
    config.bins_path       = str(tmp_path / "input_files")
    config.bin_filenames   = os.listdir(config.bins_path)
    config.reconstructions = str(tmp_path / "reconstructions")
    config.genres          = str(tmp_path / "reconstructions" / "GENREs")
    os.makedirs(config.genres)

    builder = GEMSReconstruction(config)
    builder.rast_annotate_genomes()
    builder.modelseed_reconstructions()

    assert (Path(config.genres) / "bin_6.xml").exists()
