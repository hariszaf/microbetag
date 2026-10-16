---
title: FAPROTAX functions
layout: default
parent: Modules
nav_order: 1
---

# FAPROTAX functions

The metabolic and other ecologically relevant functions that FAPROTAX assigns, grouped by theme. microbetag uses
FAPROTAX 1.2.10 by default; the pipeline can use the
[FAPROTAX v2.0 mapper](https://pages.uoregon.edu/slouca/LoucaLab/archive/FAPROTAX2/) instead (`faprotax_version: 2`,
see [Modules](modules.md)). The two versions share most functions: v2.0 adds 13, mostly fermentations, relationships
to oxygen and the degradation of polymers, and drops 4 (✓: the version has the function).

In MGG, the underscores of the names are shown as spaces.

## Methane, methanol and C1 compounds

| Function | 1.2.10 | v2.0 |
|:--|:--:|:--:|
| methanogenesis | ✓ | ✓ |
| acetoclastic_methanogenesis | ✓ | ✓ |
| hydrogenotrophic_methanogenesis | ✓ | ✓ |
| methanogenesis_by_CO2_reduction_with_H2 | ✓ | ✓ |
| methanogenesis_using_formate | ✓ | ✓ |
| methanogenesis_by_disproportionation_of_methyl_groups | ✓ | ✓ |
| methanogenesis_by_reduction_of_methyl_compounds_with_H2 | ✓ | ✓ |
| methanotrophy | ✓ | ✓ |
| methanol_oxidation | ✓ | ✓ |
| methylotrophy | ✓ | ✓ |
| reductive_acetogenesis | ✓ | ✓ |

## Chemoheterotrophy and fermentation

| Function | 1.2.10 | v2.0 |
|:--|:--:|:--:|
| chemoheterotrophy | ✓ | ✓ |
| aerobic_chemoheterotrophy | ✓ | ✓ |
| anaerobic_chemoheterotrophy |  | ✓ |
| fermentation | ✓ | ✓ |
| alcohol_fermentation |  | ✓ |
| lactic_acid_fermentation |  | ✓ |
| acetate_oxidation |  | ✓ |
| lactate_oxidation |  | ✓ |

## Relationship to oxygen

| Function | 1.2.10 | v2.0 |
|:--|:--:|:--:|
| obligate_anaerobes |  | ✓ |
| facultative_anaerobes |  | ✓ |

## Degradation of polymers and organic compounds

| Function | 1.2.10 | v2.0 |
|:--|:--:|:--:|
| cellulolysis | ✓ | ✓ |
| xylanolysis | ✓ | ✓ |
| chitinolysis | ✓ | ✓ |
| ligninolysis | ✓ | ✓ |
| pectin_degradation |  | ✓ |
| alginate_degradation |  | ✓ |
| gelatin_degradation |  | ✓ |
| collagen_degradation |  | ✓ |
| ureolysis | ✓ | ✓ |
| plastic_degradation | ✓ | ✓ |
| hydrocarbon_degradation | ✓ | ✓ |
| aliphatic_non_methane_hydrocarbon_degradation | ✓ | ✓ |
| aromatic_hydrocarbon_degradation | ✓ | ✓ |
| aromatic_compound_degradation | ✓ | ✓ |
| oil_bioremediation | ✓ | ✓ |

## Phototrophy

| Function | 1.2.10 | v2.0 |
|:--|:--:|:--:|
| phototrophy | ✓ | ✓ |
| photoautotrophy | ✓ | ✓ |
| photoheterotrophy | ✓ | ✓ |
| oxygenic_photoautotrophy | ✓ | ✓ |
| anoxygenic_photoautotrophy | ✓ | ✓ |
| anoxygenic_photoautotrophy_H2_oxidizing | ✓ | ✓ |
| anoxygenic_photoautotrophy_S_oxidizing | ✓ | ✓ |
| anoxygenic_photoautotrophy_Fe_oxidizing | ✓ | ✓ |
| aerobic_anoxygenic_phototrophy | ✓ | ✓ |
| photosynthetic_cyanobacteria | ✓ |  |
| nonphotosynthetic_cyanobacteria | ✓ |  |
| chloroplasts | ✓ |  |

## Nitrogen cycle

| Function | 1.2.10 | v2.0 |
|:--|:--:|:--:|
| nitrogen_fixation | ✓ | ✓ |
| nitrification | ✓ | ✓ |
| aerobic_ammonia_oxidation | ✓ | ✓ |
| aerobic_nitrite_oxidation | ✓ | ✓ |
| anammox | ✓ | ✓ |
| denitrification | ✓ |  |
| nitrate_denitrification | ✓ | ✓ |
| nitrite_denitrification | ✓ | ✓ |
| nitrous_oxide_denitrification | ✓ | ✓ |
| nitrate_reduction | ✓ | ✓ |
| nitrate_respiration | ✓ | ✓ |
| nitrite_respiration | ✓ | ✓ |
| nitrogen_respiration | ✓ | ✓ |
| nitrate_ammonification | ✓ | ✓ |
| nitrite_ammonification | ✓ | ✓ |

## Sulfur cycle

| Function | 1.2.10 | v2.0 |
|:--|:--:|:--:|
| respiration_of_sulfur_compounds | ✓ | ✓ |
| sulfate_respiration | ✓ | ✓ |
| sulfite_respiration | ✓ | ✓ |
| sulfur_respiration | ✓ | ✓ |
| thiosulfate_respiration | ✓ | ✓ |
| sulfur_reduction |  | ✓ |
| dark_oxidation_of_sulfur_compounds | ✓ | ✓ |
| dark_sulfide_oxidation | ✓ | ✓ |
| dark_sulfur_oxidation | ✓ | ✓ |
| dark_sulfite_oxidation | ✓ | ✓ |
| dark_thiosulfate_oxidation | ✓ | ✓ |

## Hydrogen, iron, manganese and other electron donors and acceptors

| Function | 1.2.10 | v2.0 |
|:--|:--:|:--:|
| dark_hydrogen_oxidation | ✓ | ✓ |
| knallgas_bacteria | ✓ | ✓ |
| iron_oxidation |  | ✓ |
| dark_iron_oxidation | ✓ | ✓ |
| iron_respiration | ✓ | ✓ |
| manganese_oxidation | ✓ | ✓ |
| manganese_respiration | ✓ | ✓ |
| fumarate_respiration | ✓ | ✓ |
| chlorate_reducers | ✓ | ✓ |

## Arsenic

| Function | 1.2.10 | v2.0 |
|:--|:--:|:--:|
| arsenate_detoxification | ✓ | ✓ |
| arsenate_respiration | ✓ | ✓ |
| dissimilatory_arsenate_reduction | ✓ | ✓ |
| arsenite_oxidation_detoxification | ✓ | ✓ |
| arsenite_oxidation_energy_yielding | ✓ | ✓ |
| dissimilatory_arsenite_oxidation | ✓ | ✓ |

## Hosts, symbionts and pathogens

| Function | 1.2.10 | v2.0 |
|:--|:--:|:--:|
| human_associated | ✓ | ✓ |
| human_gut | ✓ | ✓ |
| mammal_gut | ✓ | ✓ |
| human_pathogens_all | ✓ | ✓ |
| human_pathogens_diarrhea | ✓ | ✓ |
| human_pathogens_gastroenteritis | ✓ | ✓ |
| human_pathogens_meningitis | ✓ | ✓ |
| human_pathogens_nosocomia | ✓ | ✓ |
| human_pathogens_pneumonia | ✓ | ✓ |
| human_pathogens_septicemia | ✓ | ✓ |
| animal_parasites_or_symbionts | ✓ | ✓ |
| invertebrate_parasites | ✓ | ✓ |
| fish_parasites | ✓ | ✓ |
| intracellular_parasites | ✓ | ✓ |
| predatory_or_exoparasitic | ✓ | ✓ |
| plant_pathogen | ✓ | ✓ |

Of the four functions of 1.2.10 that v2.0 does not have, `chloroplasts`, `photosynthetic_cyanobacteria` and
`nonphotosynthetic_cyanobacteria` are groups of taxa rather than functions, and `denitrification` is covered by its
steps (`nitrate_denitrification`, `nitrite_denitrification`, `nitrous_oxide_denitrification`).
