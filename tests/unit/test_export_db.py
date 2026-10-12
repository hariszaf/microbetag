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

    pc = root / "pathway_compl"
    pc.mkdir()
    complement = ["md:M00019", ["K00826"], ["K01652", "K00826"], "https://kegg/url"]
    (pc / "pcompls.json").write_text(json.dumps({
        "bin_1": {"bin_1": [], "bin_2": [complement]},
        "bin_2": {"bin_1": [complement, ["md:M00020", ["K1"], ["K1", "K2"], "N/A"]], "bin_2": []},
    }))

    sc = root / "seed_compl"
    sc.mkdir()
    (sc / "seed_scores.tsv").write_text(
        "nodeA\tnodeB\tCompetitionScore\tCooperationScore\nbin_1\tbin_2\t0.59\t0.38\n")
    (sc / "all_compls.json").write_text(json.dumps({"bin_1": {"bin_2": ["cpd1", "cpd2"]}}))
    (sc / "nonseeds.json").write_text(json.dumps({"bin_1": ["cpd3"], "bin_2": ["cpd4", "cpd5"]}))
    (sc / "confidence.json").write_text(json.dumps({"bin_1": {"cpd1": 1.0}}))

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
        "genomeTaxonomy": 2, "phenDB": 2, "uniqueComplements": 2, "pathwayComplementarity": 2,
        "seedScores": 1, "SeedComplements": 1, "nonSeedSets": 2, "patricId2genomeId": 2,
    }

    # phenDB: genome ids without the FASTA extension; 33 traits, missing ones as N/A without a score (NULL)
    phen = rows(out / "phen_predictions.tsv")
    assert [r[0] for r in phen] == ["bin_1", "bin_2"]
    assert len(phen[0]) == 1 + 2 * len(export_db.PHEN_TRAITS)
    aerobe = 1 + 2 * [c for _, c in export_db.PHEN_TRAITS].index("aerobe")
    assert phen[0][aerobe:aerobe + 2] == ["YES", "0.8"]
    assert phen[1][aerobe:aerobe + 2] == ["N/A", "0.55"]
    assert phen[0][1:3] == ["N/A", "\\N"]

    # The same complement in two pairs is stored once
    assert rows(out / "unique_complements.tsv") == [
        ["1", "M00019", "K00826", "K01652;K00826"], ["2", "M00020", "K1", "K1;K2"]]
    assert rows(out / "pathway_complementarity.tsv") == [["bin_1", "bin_2", "1"], ["bin_2", "bin_1", "1,2"]]

    assert rows(out / "seed_scores.tsv") == [["bin_1", "bin_2", "0.59", "0.38"]]
    assert rows(out / "seed_complements.tsv") == [["bin_1", "bin_2", "cpd1,cpd2"]]
    assert rows(out / "genome_ids.tsv") == [["bin_1", "bin_1"], ["bin_2", "bin_2"]]
    with gzip.open(out / "all_conf.json.gz", "rt") as f:
        assert json.load(f) == {"bin_1": {"cpd1": 1.0}}


def test_cli(tmp_path, precalc, capsys):
    root, taxonomy = precalc
    export_db.main(["--precalc", str(root), "--genome-taxonomy", str(taxonomy), "--outdir", str(tmp_path / "db")])
    assert "pathwayComplementarity" in capsys.readouterr().out
