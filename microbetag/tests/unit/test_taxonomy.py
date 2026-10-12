"""Matching taxonomies to the closest genomes of a catalogue (microbetag.taxonomy)."""

import pytest

from microbetag.taxonomy import GenomeTaxonomyIndex, load_genome_taxonomy, parse_lineage

LACTIS = "d__Bacteria;p__Bacillota;c__Bacilli;o__Lactobacillales;f__Streptococcaceae;g__Lactococcus;s__Lactococcus lactis"
CATALOGUE = {
    "bin_1": LACTIS,
    "bin_2": "d__Bacteria;p__Bacillota;c__Bacilli;o__Lactobacillales;f__Streptococcaceae;g__Streptococcus;s__Streptococcus thermophilus",
    "bin_3": "d__Bacteria;p__Bacillota;c__Bacilli;o__Lactobacillales;f__Streptococcaceae;g__Streptococcus;s__Streptococcus salivarius",
    "bin_4": "d__Bacteria;p__Actinomycetota;c__Actinomycetes;o__Bifidobacteriales;f__Bifidobacteriaceae;g__Bifidobacterium;s__",
}


@pytest.fixture
def index():
    return GenomeTaxonomyIndex(CATALOGUE)


@pytest.mark.parametrize("lineage, expected", [
    (LACTIS, ("bacteria", "bacillota", "bacilli", "lactobacillales", "streptococcaceae", "lactococcus", "lactococcus lactis")),
    ("Bacteria; Bacillota ;Bacilli", ("bacteria", "bacillota", "bacilli")),
    ("D_0__Bacteria;D_1__Firmicutes", ("bacteria", "firmicutes")),
    ("Root;Bacteria;Firmicutes_A;Clostridia;unclassified_Clostridia", ("bacteria", "firmicutes_a", "clostridia")),
    ("d__Bacteria|p__Bacillota|c__|o__", ("bacteria", "bacillota")),
    ("", ()),
    (None, ()),
])
def test_parse_lineage(lineage, expected):
    assert parse_lineage(lineage) == expected


def test_same_species(index):
    m = index.closest(LACTIS)
    assert (m.genomes, m.rank) == (["bin_1"], "species")


def test_closest_is_the_genus(index):
    # A Streptococcus species that is not in the catalogue: both Streptococcus genomes
    m = index.closest("Bacteria;Bacillota;Bacilli;Lactobacillales;Streptococcaceae;Streptococcus;Streptococcus mutans")
    assert (m.genomes, m.rank) == (["bin_2", "bin_3"], "genus")


def test_formats_and_case(index):
    # SILVA-style prefixes, other case, IDTAXA's Root and unclassified ranks
    assert index.closest("D_0__bacteria;D_1__BACILLOTA;D_2__Bacilli;D_3__Lactobacillales;D_4__Streptococcaceae;"
                         "D_5__Lactococcus;D_6__Lactococcus lactis").genomes == ["bin_1"]
    m = index.closest("Root;Bacteria;Actinomycetota;Actinomycetes;Bifidobacteriales;Bifidobacteriaceae;"
                      "Bifidobacterium;unclassified_Bifidobacterium")
    assert (m.genomes, m.rank) == (["bin_4"], "genus")


def test_too_far(index):
    # Shares only the family with the catalogue: below the default minimum rank (genus)
    assert index.closest("Bacteria;Bacillota;Bacilli;Lactobacillales;Streptococcaceae;Floricoccus").genomes == []
    # ... unless the minimum rank allows it
    m = GenomeTaxonomyIndex(CATALOGUE, min_rank="family").closest(
        "Bacteria;Bacillota;Bacilli;Lactobacillales;Streptococcaceae;Floricoccus")
    assert (len(m.genomes), m.rank) == (3, "family")


def test_bare_names(index):
    assert index.closest("Lactococcus lactis").genomes == ["bin_1"]
    m = index.closest("g__Streptococcus")
    assert (m.genomes, m.rank) == (["bin_2", "bin_3"], "genus")


def test_genome_ids_are_matched_directly(index):
    assert index.match("bin_3", "anything").genomes == ["bin_3"]
    assert index.match("ASV_1", LACTIS).genomes == ["bin_1"]


def test_match_table(index):
    table = index.match_table([("bin_2", None), ("ASV_7", "Bacteria;Bacillota;Bacilli;Lactobacillales;"
                                                          "Streptococcaceae;Streptococcus")])
    assert table["genomes"].tolist() == [["bin_2"], ["bin_2", "bin_3"]]
    assert table["rank"].tolist() == ["species", "genus"]


@pytest.mark.parametrize("content", [
    "user_genome\tclassification\tfastani_reference\nbin_1\t" + LACTIS + "\tGCF_1\n",  # GTDB-Tk summary
    "genome\ttaxonomy\nbin_1\t" + LACTIS + "\n",                                           # with a header
    "bin_1\t" + LACTIS + "\n",                                                             # without a header
])
def test_load_genome_taxonomy(tmp_path, content):
    path = tmp_path / "taxonomy.tsv"
    path.write_text(content)
    assert load_genome_taxonomy(str(path)) == {"bin_1": LACTIS}
