"""
Seed sets of a GEM (PhyloMint/lib/BuildGraphNetX.py), on what they should be biologically:

- per compound, a compound is a seed or a non-seed, never both;
- the nutrients a model cannot do without, and cannot make from anything (iML1515, E. coli K-12: its metal ions
  chloride and molybdate), are seeds;
- the nutrients among the currency metabolites get a status: NH4+ and SO4 2- are non-seeds of iML1515 (made from its
  seeds: NH4+ from dopamine, SO4 2- from H2S), iron a seed (only released from imported siderophore complexes);
- compounds made from currency metabolites alone (ADP-ribose from NAD) are non-seeds.

scripts/check_seed_sets.py runs these checks, and others, on any model.
"""

import pytest

from microbetag.PhyloMint.lib.BuildGraphNetX import BIOMASS_SINK, CURRENCY_SOURCE, buildDG, getSeedSet

# Found with scripts/check_seed_sets.py: on the complete medium, iML1515 does not grow without them
IML1515_REQUIRED = {"ca2", "cl", "cobalt2", "k", "mg2", "mn2", "mobd", "ni2", "zn2"}


def seed_sets(path):
    _, seeds, nonseeds = getSeedSet(buildDG(str(path)), maxComponentSize=2)
    return set(seeds), set(nonseeds)


@pytest.mark.parametrize("model", ["bin_101.xml", "GPB:bin_000083.xml"])
def test_seeds_and_nonseeds_are_disjoint(data_dir, model):
    seeds, nonseeds = seed_sets(data_dir / "seed_compl" / "input_files" / model)
    assert seeds and nonseeds
    assert not seeds & nonseeds
    assert not {BIOMASS_SINK, CURRENCY_SOURCE} & (seeds | nonseeds)


def test_iml1515_required_nutrients_are_seeds(data_dir):
    seeds, nonseeds = seed_sets(data_dir / "seed_compl" / "iML1515.xml.gz")
    assert IML1515_REQUIRED <= seeds
    assert not seeds & nonseeds


def test_iml1515_nutrient_round_and_currency_products(data_dir):
    seeds, nonseeds = seed_sets(data_dir / "seed_compl" / "iML1515.xml.gz")
    assert {"nh4", "so4"} <= nonseeds
    assert {"fe2", "fe3"} <= seeds
    assert {"adprib", "camp", "nmn"} <= nonseeds


def test_modelseed_models_get_the_modelseed_currency():
    import cobra

    from microbetag.PhyloMint.lib.BuildGraphNetX import CURRENCY, CURRENCY_MODELSEED

    model = cobra.Model("m")
    mg_e, mg_c, atp, h2o, adp, pi, h = (cobra.Metabolite(i) for i in (
        "cpd00254_e0", "cpd00254_c0", "cpd00002_c0", "cpd00001_c0", "cpd00008_c0", "cpd00009_c0", "cpd00067_c0"))
    transport = cobra.Reaction("rxn08924_c0")  # Mg2+-importing ATPase
    transport.add_metabolites({mg_e: -1, atp: -1, h2o: -1, mg_c: 1, adp: 1, pi: 1, h: 1})
    model.add_reactions([transport])

    dg = buildDG(model)
    assert dg.graph["always_available"] == CURRENCY_MODELSEED - {"cpd00013", "cpd00048", "cpd10515", "cpd10516"}
    assert "cpd00002" not in dg and "cpd00254" in dg and not list(dg.predecessors("cpd00254"))
    assert CURRENCY != CURRENCY_MODELSEED
