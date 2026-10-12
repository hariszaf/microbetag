"""Trait prediction with phenotrex (mtg-phenotrex conda environment): genotype, then predictions."""

import os

import pytest

import microbetag
from microbetag.tools import phenotrex_genotype, phenotrex_predict

pytestmark = pytest.mark.requires_conda_env("mtg-phenotrex")


class Config:
    def __init__(self, bins, outdir):
        self.bins_path       = str(bins)
        self.phen_classes    = os.path.join(os.path.dirname(microbetag.__file__), "mtg_maps_models", "phenDB", "classes")
        self.outdir          = str(outdir)
        self.phen_traits_dir = str(outdir)
        self.genotypes_file  = str(outdir / "train.genotype")
        self.threads         = 2
        self.min_proba       = 0.6
        # Outside the container, phenotrex runs through `conda run -n mtg-phenotrex`
        self.cwd             = os.getcwd()


def test_genotype_and_predict(tmp_path, data_dir):
    config = Config(data_dir / "phenotrex" / "input_files", tmp_path)

    phenotrex_genotype(config=config)
    # header plus one line per genome
    assert len((tmp_path / "train.genotype").read_text().splitlines()) == 2

    phenotrex_predict(config=config)
    assert list(tmp_path.glob("*.prediction.tsv"))
