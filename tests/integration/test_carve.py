"""GEM reconstruction with CarveMe from bins, with FragGeneScan for gene prediction."""

from pathlib import Path

import pytest

from conftest import load_config
from microbetag.wrappers import build_genres

pytestmark = pytest.mark.requires_tools("carve", "diamond", "FragGeneScan")


def test_carve_fraggenescan(data_copy, gurobi_license):
    config = load_config(data_copy("carve") / "config_carve.yml")

    build_genres(config)

    models = {p.stem for p in Path(config.genres).glob("*.xml")}
    assert models == {"bin_101", "GPB:bin_000082"}
