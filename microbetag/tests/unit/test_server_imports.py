"""The names the server (server/services/web) imports from the package exist."""

import pytest

pytest.importorskip("mysql.connector")


def test_server_imports():
    from microbetag import Config, GetPhenotrexTraits, db, mtg_logger, run_microbetag  # noqa: F401
    from microbetag import ExportSeedComplementarities  # noqa: F401

    for name in ("GetPhenotrexTraits", "get_genomes_for_ncbi_tax_id", "get_ncbi_tax_id_for_genome",
                 "patric_from_gc_list", "get_path_compls_for_ncbi_ids", "get_path_compls_otf", "write_seed_files",
                 "seed_sets_of", "alt_genome_prefix", "build_kegg_urls"):
        assert hasattr(db, name), name


def test_server_uses_of_utils():
    """Functions of microbetag.utils the server's app_taxonomy.py uses."""
    from microbetag import utils

    assert utils.remove_nan_from_list(["GCF_1.1", float("nan"), "nan", "GCA_2.1"]) == ["GCF_1.1", "GCA_2.1"]
    assert hasattr(utils, "flatten_list")
