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
