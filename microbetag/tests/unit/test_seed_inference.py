"""
Seed sets of a GEM (PhyloMint/lib/BuildGraphNetX.py), on what they should be biologically:

- per compound, a compound is a seed or a non-seed, never both;
- the nutrients a model cannot do without, and cannot make from anything (iML1515, E. coli K-12: its metal ions
  chloride and molybdate), are seeds.

scripts/check_seed_sets.py runs these checks, and others, on any model.
"""

import pytest

from microbetag.PhyloMint.lib.BuildGraphNetX import BIOMASS_SINK, buildDG, getSeedSet

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
    assert BIOMASS_SINK not in seeds | nonseeds


def test_iml1515_required_nutrients_are_seeds(data_dir):
    seeds, nonseeds = seed_sets(data_dir / "seed_compl" / "iML1515.xml.gz")
    assert IML1515_REQUIRED <= seeds
    assert not seeds & nonseeds
