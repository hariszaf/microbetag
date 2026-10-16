"""
Pathway complementarities computed by the server from the per-genome tables (db.get_path_compls_for_ncbi_ids):
the alternatives of the beneficiary's incomplete modules whose missing KOs the donor has.
"""

import pytest

pytest.importorskip("mysql.connector")

from microbetag import db  # noqa: E402
from microbetag.pathway_complementarity import pair_complements  # noqa: E402

# uniqueComplements: id, module, missing KOs, KOs of the path
UNIQUE = [
    (1, "M00019", "K00826", "K01652;K01653;K00053;K01687;K00826"),
    (2, "M00020", "K00058;K00831", "K00058;K00831;K01079"),
]
ALTERNATIVES = {"GCF_000001.1": "1,2", "GCA_000002.1": "1"}
MODULE_KOS = [
    ("GCA_000002.1", "M00019", "K00826;K01652"),
    ("GCA_000002.1", "M00020", "K00058"),
    ("GCF_000001.1", "M00019", "K01652"),
]


def fake_select(query, values, chunk=2000):
    if "FROM genomeAlternatives" in query:
        return [(g, ids) for g, ids in ALTERNATIVES.items() if g in values]
    if "FROM genomeModuleKos" in query:
        return [row for row in MODULE_KOS if row[0] in values]
    return [row for row in UNIQUE if row[0] in values]


def test_pair_complements():
    complements = {1: {"module": "M00019", "missing": {"K00826"}},
                   2: {"module": "M00020", "missing": {"K00058", "K00831"}}}
    donor = {"M00019": {"K00826"}, "M00020": {"K00058"}}
    assert pair_complements([1, 2], donor, complements) == [1]   # M00020 also misses K00831
    assert pair_complements([2], {"M00020": {"K00058", "K00831"}}, complements) == [2]
    assert pair_complements([1], {}, complements) == []


def test_path_compls_for_ncbi_ids(monkeypatch):
    monkeypatch.setattr(db, "_select", fake_select)
    # Genomes asked for with the other prefix are found (GCA_000001.1 is stored as GCF_000001.1)
    relative = {"100": {"GCA_000001.1"}, "200": {"GCA_000002.1"}}
    res = db.get_path_compls_for_ncbi_ids(relative, {("100", "200"), ("200", "100")})

    (row,), = res[("100", "200")][("GCA_000001.1", "GCA_000002.1")]
    module, missing, path, url = row
    assert (module, missing, path) == ("M00019", "K00826", "K01652;K01653;K00053;K01687;K00826")
    assert url.startswith("https://www.kegg.jp/kegg-bin/show_pathway?map")
    # The donor GCF_000001.1 has K01652 only: no complement for GCA_000002.1
    assert ("200", "100") not in res


def test_no_genomes():
    assert db.get_path_compls_for_ncbi_ids({}, {("1", "2")}) == {}
