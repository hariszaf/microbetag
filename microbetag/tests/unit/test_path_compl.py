"""
Pathway complementarity, step by step:

1. taxon_kos_per_module(): KEGG ORTHOLOGY terms present in each genome, per KEGG module
2. all_alternatives(): for each genome and module, the KO terms missing from each alternative path
3. all_complements(): which genomes can provide the missing KO terms of another genome

These are wrapped by export_pathway_complementarities() in the main microbetag pipeline.
"""

import json

import pytest

from microbetag.helpers import MappingPaths
from microbetag.pathway_complementarity import (
    a_modules_maps,
    all_alternatives,
    all_complements,
    taxon_kos_per_module,
)
from microbetag.utils import load_merged_ko_file

N_BINS = 7


@pytest.fixture(scope="module")
def maps():
    return MappingPaths()


@pytest.fixture(scope="module")
def bin_kos_per_module(tmp_path_factory, data_dir, maps):
    # load_merged_ko_file() decompresses next to the input, so work on a copy
    import shutil
    ko_merged = tmp_path_factory.mktemp("path_compl") / "ko_merged_7bins.txt.gz"
    shutil.copy(data_dir / "path_compl" / "input_files" / "ko_merged_7bins.txt.gz", ko_merged)

    return taxon_kos_per_module(
        bins_kos_df       = load_merged_ko_file(str(ko_merged)),
        ref_ko_per_module = maps.ref_ko_per_module,
    )


def test_taxon_kos_per_module(bin_kos_per_module):
    assert len(bin_kos_per_module) == N_BINS


def test_all_alternatives_and_complements(tmp_path, maps, bin_kos_per_module):
    alts_file = tmp_path / "alternatives.json"
    pc_file   = tmp_path / "complementarities.json"

    bins_alternatives = all_alternatives(
        bin_kos_per_module           = bin_kos_per_module,
        modules_definitions_json_map = maps.modules_definitions_json_map,
        alts_output_file             = str(alts_file),
    )
    assert alts_file.exists()
    assert set(bins_alternatives) <= set(bin_kos_per_module)

    complements = all_complements(
        bin_kos_per_module,
        bins_alternatives,
        a_modules_maps(maps.kegg_modules_to_maps),
        str(pc_file),
        tinyurl=False,
    )
    assert set(complements) == set(bin_kos_per_module)

    with open(pc_file) as f:
        assert len(json.load(f)) == N_BINS
