"""
Build the microbetag-annotated CX2 network from the outputs of the other steps
(FAPROTAX, phenotrex, pathway and seed complementarities, manta) stored in tests/data/build_cx2.
"""

import pytest

from conftest import KOFAM_DIR, load_config
from microbetag.build_mtg_cx2 import mtg_annotate_network


@pytest.mark.xfail(
    reason="tests/data/build_cx2/output_files predates the current seed complementarity outputs "
           "(no kegg_module_related_*.pkl); regenerate them with a full run",
    raises=FileNotFoundError,
    strict=True,
)
def test_build_cx2(data_copy):
    out = data_copy("build_cx2")
    before = set(out.rglob("*.cx2"))

    # Pathway complementarities are precomputed, but the config still needs a kofam folder
    config = load_config(out / "config_buildCX2.yml", kofam_database={"dir_path": str(KOFAM_DIR)})
    mtg_annotate_network(config)

    assert set(out.rglob("*.cx2")) - before
