"""
Aim:
    This script gets the KEGG compounds for each KEGG module
    and then maps them to their corresponding ModelSEED compound ids.
    Last, it writes a 3-column file called `` whith the KEGG module id, the 
    KEGG compound id and its corresponding ModelSEED id.

Author:
    Haris Zafeiropoulos

Notes:

    Run using the modelseed conda env in Genius

    For the following KEGG modules, no KEGG compound involved was mapped to a ModelSEED compound
    "M00065",  # Glycan biosynthesis:  GPI-anchor biosynthesis
    "M00068",  # Glycan biosynthesis:  Glycosphingolipid biosynthesis: globo-series
    "M00069",  # Glycan biosynthesis:  Glycosphingolipid biosynthesis: ganglio series
    "M00070",  # Glycan biosynthesis:  Glycosphingolipid biosynthesis: lacto-series
    "M00071",  # Glycan biosynthesis:  Glycosphingolipid biosynthesis: neolacto-series
    "M00831",  # Biosynthesis of terpenoids and polyketides:  Kedarcidin 2-hydroxynaphthoate moiety biosynthesis
    "M00832",  # Biosynthesis of terpenoids and polyketides:  Kedarcidin 2-aza-3-chloro-beta-tyrosine moiety biosynthesis
    "M00867"   # Lipopolysaccharide metabolism:  KDO2-lipid A modification pathway

    mean_non_seedset_length: 752.8452673677974
    mean_seedset_length: 341.1272404088283
"""


"""
Part A: 
    - remove compounds from seed sets that are related to environmental metabolites that can be produced in several ways.
    - remove from non seed sets compounds that cannot be produced in any other way than from entering the cell from the environment.
"""
import os
import sys
import time
import cobra
import json


# for i in range(0,12,2):
#     subfolders = [ "subfolder_" + str(i) for i in range(i,i+2) ]
#     print(subfolders)

# subfolders = [ "subfolder_" + str(i) for i in range(8,12) ]

f = open("stats2.tsv", "w")
f.write("model_id" + "\t" + "environmental init seeds" + "\t" + "non environmental init seeds" +  "\t" +  "total init seeds" + "\t" + "updated seeds" + "\t" +  "initial non seeds" + "\t" + "updated non seeds" + "\n")
pwd = os.getcwd()
xmls_path = "/home/luna.kuleuven.be/u0156635/Documents/projects/microbetag/gtdb_modelseed_gems/"
dics_path = "/home/luna.kuleuven.be/u0156635/Documents/projects/microbetag/dics/"

# subfolders = ['subfolder_0', 'subfolder_1']
# subfolders = ['subfolder_2', 'subfolder_3']
subfolders = ['subfolder_4', 'subfolder_5']
# subfolders = ['subfolder_6', 'subfolder_7'] 
# subfolders = ['subfolder_8', 'subfolder_9']
subfolders = ['subfolder_10', 'subfolder_11']

parsed_ids = set()
for subfolder in subfolders:
    updated_seeds = {}
    updated_nonSeeds = {}
    current_seeds = json.load(open(dics_path + subfolder + "_ConfDic.json", "r"))
    current_nonSeeds = json.load(open(dics_path + subfolder + "_nonSeedSetDic.json", "r"))
    for xml in os.listdir("/".join([xmls_path, subfolder])):
        print("xml:", xml)
        s1 = time.time()
        counter = 0; counter2 = 0
        xml_path = "/".join([xmls_path, subfolder, xml])
        command = " ".join(["gunzip", xml_path])
        model_file = xml_path[:-3]
        model_id = xml[:-7]
        if model_id in parsed_ids:
            print("model already parsed")
            continue
        os.system(command)
        model = cobra.io.read_sbml_model(model_file)
        model_mets = [met.id for met in model.metabolites]
        models_tmp_non_seeds = current_nonSeeds[model_id]
        models_tmp_seeds = current_seeds[model_id]
        models_seeds = []; models_nonSeeds = models_tmp_non_seeds.copy()
        for pot_seed in models_tmp_seeds:
            check = True
            if "_e0" in pot_seed:
                counter += 1
                cor_in_met = pot_seed[:-3] + "_c0"
                if cor_in_met in models_tmp_seeds:
                     # Both _c0 and _e0 among the potential seed set.
                     models_seeds.append(pot_seed)
                else:
                    cor_in_met = cor_in_met[2:]
                    if cor_in_met not in model_mets:
                        # The _c0 case is not among the model's metabolites.
                        models_seeds.append(pot_seed)
                    else:
                        for rxn in model.metabolites.get_by_id(cor_in_met).reactions:
                            if cor_in_met in [met.id for met in rxn.products]:
                                if pot_seed[2:] not in [met.id for met in rxn.reactants]:
                                    # There is at least a reaction that does not include the _e0 case that produces the _c0 metabolite.
                                    check = False
                                    break
                        if check:
                            models_seeds.append(pot_seed)
                            models_nonSeeds.remove("M_"+cor_in_met)
            else:
                counter2 += 1
                models_seeds.append(pot_seed)
        with open("stats3.tsv", "a") as f:
            f.write(model_id + "\t" + str(counter) + "\t" + str(counter2) + "\t" + str(len(models_tmp_seeds)) + "\t" + str(len(models_seeds)) + "\t" + str(len(models_tmp_non_seeds)) + "\t" + str(len(models_nonSeeds)) + "\n")
        updated_seeds[model_id] = models_seeds
        updated_nonSeeds[model_id] = models_nonSeeds
        command = " ".join(["gzip", model_file])
        os.system(command)
        s2 = time.time()
        parsed_ids.add(model_id)
        print(str(s2-s1), "seconds for a .xml")
    with open("updated_" + subfolder + "_SeedsDic.json", "w") as f:
        json.dump(updated_seeds, f)
    with open("updated_" + subfolder + "_NonSeedsDic.json", "w") as f:
        json.dump(updated_nonSeeds, f)




"""
Part B:
    - Map KEGG terms to their corresponding ModelSEED compounds.
"""
modelseed_compounds_metadata = open("modelseedDB_compounds.tsv", "r")
modelseedId_to_keggId = {}
keggId_to_modelSeedId = {}
counter = 0
all_lines = 0
for line in modelseed_compounds_metadata:
    all_lines += 1
    modelseedId = line.split("\t")[0]
    keggId_parts = line.split("\t")[18].split("KEGG")
    if len(keggId_parts) > 1:
        keggId = keggId_parts[1][2:]
        if "|" in keggId:
            keggId = keggId.split("|")[0]
        if ";" in keggId:
            modelseedId_to_keggId[modelseedId] = keggId
            keggIds = keggId.split(";")
            for keggId in keggIds:
                if keggId not in keggId_to_modelSeedId:
                    keggId_to_modelSeedId[keggId] = modelseedId
                else:
                    keggId_to_modelSeedId[keggId] += ";" + modelseedId
        else:
            if keggId not in keggId_to_modelSeedId:
                keggId_to_modelSeedId[keggId] = modelseedId
            else:
                keggId_to_modelSeedId[keggId] += ";" + modelseedId
            modelseedId_to_keggId[modelseedId] = keggId
    else:
        counter += 1

print("ModelSEED compounds with not a corresponding KEGG one:", str(counter))
print("ModelSEED compound ids in total:", all_lines)


"""
Part C:
    - For the KEGG modules of interest, get the KEGG terms that participate in it and map them to their corresponding ModelSEED ids
    ModelSEED   KEGG COMPOUND   KEGG MODULE
    cpd00020        C00022      M00001 
    cpd00061        C00074      M00001 
"""
import requests

modules_file = open("module_map_pairs.tsv", "r")
modules = [x.split("\t")[0][3:] for x in modules_file.readlines() ]
base_url = "https://rest.kegg.jp/link/cpd/"
f = open("seedId_keggId_module.tsv", "w")
counter2 = 0
for module in modules:
    print("Processing module...", module)
    counter2 += 1
    if counter2 % 10 == 0:
        print(counter2, " out of ", str(len(modules)), " modules.")
    time.sleep(1)
    url = "".join([base_url, module])
    r = requests.get(url)
    modules_kegg_compounds = r.text.split("\n")
    for line in modules_kegg_compounds:
        try:
            kegg_compound = line.split("\t")[1][4:]
        except:
            pass
        print(kegg_compound)
        if kegg_compound in keggId_to_modelSeedId:
            modelseed_compound = keggId_to_modelSeedId[kegg_compound]
            if ";" in modelseed_compound:
                modelseed_compounds = modelseed_compound.split(";")
                for modelseed_compound in modelseed_compounds:
                    f.write("\t".join( [modelseed_compound, kegg_compound, module]))
                    f.write("\n")
            else:
                f.write("\t".join( [modelseed_compound, kegg_compound, module]))
                f.write("\n")


"""
Part D: 
    In the seedId_keggId_module.tsv file we have saved the compounds of interest both in ModelSEED and in KEGG terms
    Here, we use this file and we find the overlaps between each species seed and non seed sets with the compounds of interest. 

"""
import json
import pandas as pd
import pickle
import os

modules_compounds = pd.read_csv("seedId_keggId_module.tsv", sep="\t")
modules_compounds.columns = ["modelseed", "kegg", "module"]
modelseed_compounds_of_interest = set(modules_compounds["modelseed"].unique().tolist())

print(modelseed_compounds_of_interest)
print(len(modelseed_compounds_of_interest))

number_of_models = 0
patricId_to_seeds_of_interest = {}
patricId_to_non_seeds_of_interest = {}

mean_non_seedset_length = 0 ; mean_non_seedset_of_interest = 0
non_seedset_files = [file for file in os.listdir(".") if file.endswith("_NonSeedsDic.json")]  ## nonSeedSetDic
for sfile in non_seedset_files:
    non_seedset_file = json.load(open(sfile,"r"))
    model_names = list(non_seedset_file.keys())
    for smodel_name in model_names:
        non_seedset = set( [x[2:] for x in non_seedset_file[smodel_name]] )
        non_seedset_no_compartments = set( [x[2:-3]  for x in non_seedset_file[smodel_name] ])
        non_seeds_of_interest = non_seedset_no_compartments.intersection(modelseed_compounds_of_interest)
        mean_non_seedset_length += len(non_seedset)  # stat
        non_seeds_of_interest_tmp = list(non_seeds_of_interest.copy())
        patricId_to_non_seeds_of_interest[smodel_name] = non_seeds_of_interest
        mean_non_seedset_of_interest += len(non_seeds_of_interest)  # stat
        """
        for y in non_seeds_of_interest:
            if "".join([y, "_c0"]) not in non_seedset:
                index = non_seeds_of_interest_tmp.index(y)
                del non_seeds_of_interest_tmp[index]
        patricId_to_non_seeds_of_interest[model_name] = set(non_seeds_of_interest_tmp)
        mean_non_seedset_of_interest += len(non_seeds_of_interest_tmp)
        """

seedset_files = [file for file in os.listdir(".") if file.endswith("_SeedsDic.json")]  # ConfDic.json
mean_seedset_length = 0
mean_seedset_of_interest = 0

for sfile in seedset_files:
    seedset_file = json.load(open(sfile,"r"))
    model_names = list(seedset_file.keys())
    for model_name in model_names:
        number_of_models += 1
        # Get seeds with and without their compartment specific part
        seedset = set([x[2:]  for x in seedset_file[model_name]])  # x[2:-3] seedset_file[model_name].keys()
        seedset_no_compartments = set([x[2:-3]  for x in seedset_file[model_name]])  # seedset_file[model_name].keys()
        mean_seedset_length += len(seedset)
        # Using the list with the compounds without their compartment, find which are potentlially of interest
        seeds_of_interest = seedset_no_compartments.intersection(modelseed_compounds_of_interest)
        seeds_of_interest_tmp = list(seeds_of_interest.copy())
        for pot_seed in seeds_of_interest:
            if pot_seed in patricId_to_non_seeds_of_interest[model_name]:
                seeds_of_interest_tmp.remove(pot_seed)
        patricId_to_seeds_of_interest[model_name] = set(seeds_of_interest_tmp)
        mean_seedset_of_interest += len(set(seeds_of_interest_tmp))
        """
        seeds_of_interest_tmp = list(seeds_of_interest.copy())
        # Check whether the potential seed is also part of the cytosol and not just an exchange reaction
        for x in seeds_of_interest:
            # print(x)
            # print("".join([x, "_c0"]))
            if "".join([x, "_c0"]) not in seedset:
                index = seeds_of_interest_tmp.index(x)
                del seeds_of_interest_tmp[index]
        # Assign the list of seeds to a model
        patricId_to_seeds_of_interest[model_name] = set(seeds_of_interest_tmp)
        mean_seedset_of_interest += len(seeds_of_interest_tmp)
        """


print("mean length of initial seedset:", str(mean_seedset_length/number_of_models))
print("mean length of seedsets of interest:", str(mean_seedset_of_interest/number_of_models))
print("~~")
print("mean of initial length of non seed sets:", str(mean_non_seedset_length/number_of_models))
print("mean of non seed sets of interest:", str(mean_non_seedset_of_interest/number_of_models))

tmp_dict = {key: list(value) for key, value in patricId_to_seeds_of_interest.items()}
df1 = pd.DataFrame(list(tmp_dict.items()), columns=['PATRIC', 'SeedSet'])
df1['PATRIC'] = df1['PATRIC'].str.replace('.PATRIC', '')
df1.set_index('PATRIC', inplace=True)
with open("seedsets_of_interest.pckl","wb") as f:
    # tmp_dict = {key: list(value) for key, value in patricId_to_seeds_of_interest.items()}
    # df1 = pd.DataFrame(list(tmp_dict.items()), columns=['PATRIC', 'SeedSet'])
    pickle.dump(df1, f)

tmp_dict = {key: list(value) for key, value in patricId_to_non_seeds_of_interest.items()}
df2 = pd.DataFrame(list(tmp_dict.items()), columns=['PATRIC', 'NonSeedSet'])
df2['PATRIC'] = df2['PATRIC'].str.replace('.PATRIC', '')
df2.set_index('PATRIC', inplace=True)
with open("non_seedsets_of_interest.pckl","wb") as f:
    # tmp_dict = {key: list(value) for key, value in patricId_to_non_seeds_of_interest.items()}
    # df2 = pd.DataFrame(list(tmp_dict.items()), columns=['PATRIC', 'NonSeedSet'])
    pickle.dump(df2, f)


"""
PART E:
    for each species and their seed sets
    get their oversection with each other species' nonSeedSets
    (thus, the latter will provide a seed node to the beneficiary )
"""







