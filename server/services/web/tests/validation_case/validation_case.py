import time

import pandas as pd
from microbetag import microbetag

variovorax_ids = pd.read_csv(
    "tests/validation_case/data/variovorax_ncbiIds.csv", header=None
)
arnimatodetes_ids = pd.read_csv(
    "tests/validation_case/data/arnimatodetes_ncbiIds.csv", header=None
)
chloroflexi_ids = pd.read_csv(
    "tests/validation_case/data/chloroflexi_ncbiIds.csv", header=None
)
kapabacteris_ids = pd.read_csv(
    "tests/validation_case/data/kapabacteria_ncbiIds.csv", header=None
)
saccharibacteria_ids = pd.read_csv(
    "tests/validation_case/data/saccharibacteria_ncbiIds.csv", header=None
)
microbacterium_ids = pd.read_csv(
    "tests/validation_case/data/microbacterium_ncbiIds.csv", header=None
)

variovorax_ids.columns = (
    chloroflexi_ids.columns
) = (
    arnimatodetes_ids.columns
) = (
    kapabacteris_ids.columns
) = saccharibacteria_ids.columns = microbacterium_ids.columns = ["ncbiId"]

nonseeds, kmap = microbetag.load_seed_complement_files()


def get_variable_name(value, namespace):
    for name, obj in namespace.items():
        if obj is value:
            return name
    return None


def get_compl_between_vario_and(interaction_species_ids):
    """
    Exports complements between Variovorax and a taxonomic group.
    Both correspond to a list of NCBI Taxonomy ids.
    The function returns 2 dictionaries with keys like:
        ncbiIdForVorovoraxSpecies_ncbiIdForTaxonSpecies
    and the other way around and values a list with the complements.
    """
    namespace = locals()
    interaction_taxon_name = get_variable_name(interaction_species_ids, namespace)
    vario_benefits_from = {}
    vario_provides_to = {}
    counter = 0
    for index_beneficiary, varioId in variovorax_ids.iterrows():
        t1 = time.time()
        counter += 1
        # print(counter)
        # if counter <4:
        #     continue
        vvarioId = str(varioId["ncbiId"])
        print(
            "\n~~~~~\n>>>>>",
            str(counter),
            " out of ",
            str(len(variovorax_ids)),
            " Variovorax NCBI ids.",
        )
        inner_counter = 0
        for index_donor, arnimaId in interaction_species_ids.iterrows():
            inner_counter += 1
            print(
                str(inner_counter),
                "out of",
                str(len(interaction_species_ids)),
                interaction_taxon_name,
                "NCBI ids.",
            )
            aarnimaId = str(arnimaId["ncbiId"])
            try:
                seed_compl = microbetag.get_paired_seed_complements(
                    beneficiary=vvarioId,
                    donor=aarnimaId,
                    kmap=kmap,
                    nonseeds=nonseeds,
                    type="ncbiTaxonomyIds",
                )
            except ValueError:
                print("No good", vvarioId, aarnimaId)
                pass
            akey = "_".join([vvarioId, aarnimaId])
            vario_benefits_from[akey] = seed_compl
            try:
                seed_compl = microbetag.get_paired_seed_complements(
                    beneficiary=aarnimaId,
                    donor=vvarioId,
                    kmap=kmap,
                    nonseeds=nonseeds,
                    type="ncbiTaxonomyIds",
                )
                akey = "_".join([aarnimaId, vvarioId])
                vario_provides_to[akey] = seed_compl
            except ValueError:
                print("No good", aarnimaId, vvarioId)
                pass
        t2 = time.time()
        print("For an outer loop: ", str(t2 - t1), "seconds.")
    return vario_benefits_from, vario_provides_to


def ratio_of_pairs_of_models_checked_vs_theoritical(
    ncbiListA,
    ncbiListB,
):
    """
    If no GEMs present for a pair of NCBI ids we get a 'no pair of GEMs' message.

    NOTE: not sure why we need this... ? the return does not make sense for the time being
    """
    theoretical = len(ncbiListA) * len(ncbiListB)

    return theoretical


def get_unique_maps_along_with_their_complements(dict_to_check):
    unique_maps = {}
    for ncbi_pair, set_of_compl in dict_to_check.items():
        for compl in set_of_compl:
            if isinstance(compl, list):
                if compl[1] not in unique_maps:
                    unique_maps[compl[1]] = set()
                    unique_maps[compl[1]].add(compl[-1])
                else:
                    unique_maps[compl[1]].add(compl[-1])
    return unique_maps


def get_unique_maps_for_a_met_category(dict_to_check, name_of_the_map):
    """
    Get unique complements related with a specific map included in a beneficiayr-donor dictionary.
    Returns a set with URLs of colored KEGG maps.
    """
    unique_maps = set()
    for ncbi_pair, set_of_compl in dict_to_check.items():
        for compl in set_of_compl:
            if isinstance(compl, list):
                if name_of_the_map in compl[1]:
                    unique_maps.add(compl[-1])
    return unique_maps


(
    vario_benefits_from_kapabacteria,
    kapabacteria_benefit_from_vario,
) = get_compl_between_vario_and(kapabacteris_ids)
(
    vario_benefits_from_saccharibacteria,
    saccharibacteria_benefit_from_vario,
) = get_compl_between_vario_and(saccharibacteria_ids)


(
    vario_benefits_from_microbacterium,
    microbacterium_benefits_from_vario,
) = get_compl_between_vario_and(microbacterium_ids)
