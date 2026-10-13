"""
Seed complementarity with ExportSeedComplementarities, in three cases:

1. seed and non-seed sets computed from the GEMs, then scores and complements
2. scores and complements from previously computed sets (the reference ones in tests/data)
3. seed and non-seed sets only, kept in the BiGG namespace
"""

import json

import pytest

from microbetag.helpers import MappingPaths
from microbetag.seed_complementarity import ExportSeedComplementarities

# Reference scores for the two test GEMs: beneficiary, donor, competition, cooperation
# (seed sets per compound: PhyloMint/lib/BuildGraphNetX.py, see "How microbetag gets seed sets" in the README)
EXPECTED_SCORES = {
    "GPB:bin_000083\tbin_101\t0.39\t0.48",
    "bin_101\tGPB:bin_000083\t0.31\t0.27",
}
MODELS = {"GPB:bin_000083", "bin_101"}


class Config:
    def __init__(self, models_dir, out_dir):
        maps = MappingPaths()

        self.genres = self.for_reconstructions = str(models_dir)
        self.get_scores      = True
        self.get_complements = True

        self.seeds_outdir    = str(out_dir)
        self.sc_pkl          = str(out_dir / "seed_complements.pckl")
        self.module_seeds    = str(out_dir / "kegg_module_related_seeds.pckl")
        self.module_nonseeds = str(out_dir / "kegg_module_related_nonseeds.pckl")

        self.seed_ko_mo                = maps.seed_ko_mo
        self.genre_reconstruction_with = "carveme"
        self.metanetx_compounds        = maps.metanetx_compounds

        # Map BiGG compounds to ModelSEED
        self.switch_namespace = True

        self.threads       = 2
        self.skip_sets     = False
        self.prev_conf     = None
        self.prev_nonseeds = None


@pytest.fixture
def seed_data(data_dir):
    return data_dir / "seed_compl"


def read_scores(out_dir):
    return set((out_dir / "seed_scores.tsv").read_text().splitlines())


def test_sets_scores_and_complements(tmp_path, seed_data):
    ss = ExportSeedComplementarities(config=Config(seed_data / "input_files", tmp_path))
    ss.get_sets()
    ss.get_scores_and_compls()

    with open(tmp_path / "SeedSetDic.json") as f:
        assert set(json.load(f)) == MODELS
    assert read_scores(tmp_path) == EXPECTED_SCORES
    assert (tmp_path / "seed_complements.pckl").exists()


def test_scores_from_previous_sets(tmp_path, seed_data):
    config               = Config(seed_data / "input_files", tmp_path)
    config.skip_sets     = True
    config.prev_conf     = str(seed_data / "output_files" / "confidenceDic.json")
    config.prev_nonseeds = str(seed_data / "output_files" / "nonSeedSetDic.json")

    ss = ExportSeedComplementarities(config=config)
    ss.get_sets()
    ss.get_scores_and_compls()

    assert read_scores(tmp_path) == EXPECTED_SCORES


def test_bigg_sets_only(tmp_path, seed_data):
    config                  = Config(seed_data / "input_files", tmp_path)
    config.get_complements  = False
    config.switch_namespace = False

    ExportSeedComplementarities(config=config).get_sets()

    with open(tmp_path / "SeedSetDic.json") as f:
        result = json.load(f)
    with open(seed_data / "output_files_bigg" / "SeedSetDic.json") as f:
        expected = json.load(f)
    assert {k: set(v) for k, v in result.items()} == {k: set(v) for k, v in expected.items()}
