"""Exporting pre-calculations as the files of a microbetag database (microbetag-export-db)."""

import gzip
import json

import pytest

from microbetag.cli import export_db

LINEAGE = "d__Bacteria;p__Bacillota;c__Bacilli;o__Lactobacillales;f__Streptococcaceae;g__Lactococcus;s__"


@pytest.fixture
def precalc(tmp_path):
    """A small precalc output folder, as the microbetag-nf precalc workflow writes it."""
    root = tmp_path / "precalc"

    pred = root / "phenotrex" / "predictions"
    pred.mkdir(parents=True)
    (pred / "aerobe.tsv").write_text(
        "# Trait: aerobe\nIdentifier\tTrait present\tConfidence\nbin_1.fa\tYES\t0.8\nbin_2.fasta\tN/A\t0.55\n")

    # The alternatives of the incomplete modules of each genome, and its KOs (precalc's alternatives.json and
    # ko_merged.txt): the same alternative of M00019 in both genomes is stored once
    pc = root / "pathway_compl"
    pc.mkdir()
    valine = {"['K01652', 'K00826']": ["K00826"]}
    (pc / "alternatives.json").write_text(json.dumps({
        "bin_1": {"md:M00019": valine},
        "bin_2": {"md:M00019": valine, "md:M00020": {"['K00058', 'K00831']": ["K00831", "K00058"]}},
    }))
    hmm = root / "hmmsearch"
    hmm.mkdir()
    (hmm / "ko_merged.txt").write_text(
        "bin_id\tcontig_id\tko_term\nbin_1\tK00826\tc1_1\nbin_1\tK00058\tc1_2\nbin_2\tK01652\tc2_1\n")

    sc = root / "seed_compl"
    sc.mkdir()
    (sc / "nonseeds.json").write_text(json.dumps({"bin_1": ["cpd3"], "bin_2": ["cpd4", "cpd5"]}))
    (sc / "confidence.json").write_text(json.dumps({"bin_1": {"cpd1": 1.0}}))
    (sc / "keggm_nonseeds.pkl").write_bytes(b"pickle")

    taxonomy = tmp_path / "taxonomy.tsv"
    taxonomy.write_text(f"bin_1\t{LINEAGE}\nbin_2\t{LINEAGE}\n")
    return root, taxonomy


def rows(path):
    return [line.split("\t") for line in path.read_text().splitlines()]


def test_export(tmp_path, precalc):
    root, taxonomy = precalc
    out = tmp_path / "db"

    counts = export_db.export(str(root), str(taxonomy), str(out))

    assert counts == {
        "genomeTaxonomy": 2, "phenDB": 2, "uniqueComplements": 2, "genomeAlternatives": 2, "genomeModuleKos": 7,
        "seed sets": 1, "patricId2genomeId": 2, "genome2taxNcbiId": 0,
    }

    # phenDB: genome ids without the FASTA extension; 33 traits, missing ones as N/A without a score (NULL)
    phen = rows(out / "phen_predictions.tsv")
    assert [r[0] for r in phen] == ["bin_1", "bin_2"]
    assert len(phen[0]) == 1 + 2 * len(export_db.PHEN_TRAITS)
    aerobe = 1 + 2 * [c for _, c in export_db.PHEN_TRAITS].index("aerobe")
    assert phen[0][aerobe:aerobe + 2] == ["YES", "0.8"]
    assert phen[1][aerobe:aerobe + 2] == ["N/A", "0.55"]
    assert phen[0][1:3] == ["N/A", "\\N"]

    # Each alternative once; the alternatives of each genome; the KOs of each genome per module (M00019:
    # K01652, K00826; M00020: K00058, K00831; the KOs are in other modules too)
    assert rows(out / "unique_complements.tsv") == [
        ["1", "M00019", "K00826", "K01652;K00826"], ["2", "M00020", "K00058;K00831", "K00058;K00831"]]
    assert rows(out / "genome_alternatives.tsv") == [["bin_1", "1"], ["bin_2", "1,2"]]
    assert [r for r in rows(out / "genome_module_kos.tsv") if r[1] in ("M00019", "M00020")] == [
        ["bin_1", "M00019", "K00826"], ["bin_1", "M00020", "K00058"], ["bin_2", "M00019", "K01652"]]
    assert not (out / "pathway_complementarity.tsv").exists()

    assert rows(out / "genome_ids.tsv") == [["bin_1", "bin_1"], ["bin_2", "bin_2"]]
    with gzip.open(out / "all_conf.json.gz", "rt") as f:
        assert json.load(f) == {"bin_1": {"cpd1": 1.0}}
    assert (out / "all_nonseeds_module_related.pckl").read_bytes() == b"pickle"


def test_cli(tmp_path, precalc, capsys):
    root, taxonomy = precalc
    export_db.main(["--precalc", str(root), "--genome-taxonomy", str(taxonomy), "--outdir", str(tmp_path / "db")])
    assert "genomeAlternatives" in capsys.readouterr().out


def test_gem_and_ncbi_ids(tmp_path, precalc):
    root, taxonomy = precalc
    (tmp_path / "gems.tsv").write_text("bin_1\t100.11\nbin_2\t100.12\n")
    (tmp_path / "ncbi.tsv").write_text("1358\tbin_1\n")
    out = tmp_path / "db"
    export_db.export(str(root), str(taxonomy), str(out), str(tmp_path / "gems.tsv"), str(tmp_path / "ncbi.tsv"))
    assert rows(out / "genome_ids.tsv") == [["bin_1", "100.11"], ["bin_2", "100.12"]]
    assert rows(out / "genome_ncbi_tax_ids.tsv") == [["1358", "bin_1"]]
