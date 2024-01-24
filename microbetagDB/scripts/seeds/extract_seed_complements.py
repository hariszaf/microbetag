"""
Aim:

Author:
    Haris Zafeiropoulos


"""
import sys

mhelp = """
        You need to provide the script with an argument describing which part you need to implement by typing one of the following:
        update
        ofInterest
        complements     to extract complements among seed and non seed sets
        format  to
"""

if len(sys.argv) < 2:
    print(mhelp)
    sys.exit(0)


if sys.argv[1] == "update":
    """
    Part A:
        Update seed sets returned by PhyloMint by:
        - removing compounds from seed sets that are related to environmental metabolites that can be produced in several ways within the cell.
        - removing from non seed sets compounds that cannot be produced in any other way than from entering the cell from the environment.
    """
    import cobra
    import json
    import os
    import time

    f = open("stats31.tsv", "w")
    f.write("model_id" + "\t" + "environmental init seeds" + "\t" + "non environmental init seeds" +  "\t" +  "total init seeds" + "\t" + "updated seeds" + "\t" +  "initial non seeds" + "\t" + "updated non seeds" + "\n")
    pwd = os.getcwd()
    xmls_path = "/home/luna.kuleuven.be/u0156635/Documents/projects/microbetag/gtdb_modelseed_gems/"
    dics_path = "/home/luna.kuleuven.be/u0156635/Documents/projects/microbetag/dics/initial/"

    # NOTE: For computing resources related reasons, we had to iterate using up to 2 cases each time
    # subfolders = ['subfolder_0', 'subfolder_1']
    # subfolders = ['subfolder_2', 'subfolder_3']
    # subfolders = ['subfolder_4', 'subfolder_5']
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
            with open("stats31.tsv", "a") as f:
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

elif sys.argv[1] == "ofInterest":
    """
    Part B:
        In the seedId_keggId_module.tsv file we keep the compounds of interest both in ModelSEED
        and in KEGG terms, meaning compounds that participate at least in one KEGG MODULE.
        Here, using this file we keep extract those compounds of interest out of each species seed and non seed sets.
    """
    import json
    import os
    import pandas as pd
    import pickle

    base_path = "/".join(os.path.abspath(__file__).split("/")[:-4])
    modules_compounds = pd.read_csv(base_path + "/microbetagDB/mappings/kegg_mappings/seedId_keggId_module.tsv", sep="\t")
    modules_compounds.columns = ["modelseed", "kegg", "module"]
    modelseed_compounds_of_interest = set(modules_compounds["modelseed"].unique().tolist())

    print(len(modelseed_compounds_of_interest))

    number_of_models = 0
    patricId_to_seeds_of_interest = {}
    patricId_to_non_seeds_of_interest = {}

    mean_non_seedset_length = 0 ; mean_non_seedset_of_interest = 0
    dics_path  = "/home/luna.kuleuven.be/u0156635/Documents/projects/microbetag/dics/updated/"
    non_seedset_files = [file for file in os.listdir(dics_path) if file.endswith("_NonSeedsDic.json")]  ## nonSeedSetDic
    for sfile in non_seedset_files:
        non_seedset_file = json.load(open(dics_path + sfile,"r"))
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

    seedset_files = [file for file in os.listdir(dics_path) if file.endswith("_SeedsDic.json")]  # ConfDic.json
    mean_seedset_length = 0
    mean_seedset_of_interest = 0

    for sfile in seedset_files:
        seedset_file = json.load(open(dics_path + sfile,"r"))
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

    with open("updated_seedsets_of_interest2.pckl","wb") as f:
        pickle.dump(df1, f)

    tmp_dict = {key: list(value) for key, value in patricId_to_non_seeds_of_interest.items()}
    df2 = pd.DataFrame(list(tmp_dict.items()), columns=['PATRIC', 'NonSeedSet'])
    df2['PATRIC'] = df2['PATRIC'].str.replace('.PATRIC', '')
    df2.set_index('PATRIC', inplace=True)

    with open("updated_non_seedsets_of_interest2.pckl","wb") as f:
        pickle.dump(df2, f)


elif sys.argv[1] == "complements":
    """
    PART C:
        for each species and their seed sets
        get their oversection with each other species' nonSeedSets
        (thus, the latter will provide a seed node to the beneficiary )
    """
    import pandas as pd
    from joblib import Parallel, delayed
    from tqdm import tqdm
    import pickle

    with open("updated_seedsets_of_interest.pckl", "rb") as f:
        df1 = pickle.load(f)
    with open("updated_non_seedsets_of_interest.pckl", "rb") as f:
        df2 = pickle.load(f)
    with open("updated_seedsets_of_interest_no3.pckl", "rb") as f:
        df_no3 = pickle.load(f)

    seeds_3000_df = df1.merge(df_no3, on='PATRIC', how='left', indicator=True).query('_merge == "left_only"').drop('_merge', axis=1)
    seeds_3000_df = seeds_3000_df.drop("SeedSet_y", axis =1)
    seeds_3000_df.columns = ["SeedSet"]

    nonseeds_3000_df = df2.merge(df_no3, on='PATRIC', how='left', indicator=True).query('_merge == "left_only"').drop('_merge', axis=1)
    nonseeds_3000_df = nonseeds_3000_df.drop("SeedSet", axis=1)

    # NOTE: For computing resources related reasons we had to implement this task in chunks
    # df1_subset = df1.head(2000)
    # df2_subset = df1.iloc[2000:4000]
    # df3_subset = df1.iloc[4000:6000]
    # df4_subset = df1.iloc[8000:10000]
    # df5_subset = df1.iloc[10000:12000]
    # df6_subset = df1.iloc[12000:18000]
    # df9_subset = df1.iloc[18000:20000]
    # df10_subset = df1.iloc[20000:22000]
    # df11_subset = df1.iloc[22000:24000]
    # df12_subset = df1.iloc[24000:26000]
    # df13_subset = df1.iloc[26000:28000]
    df14_subset = df1.iloc[28000:]

    # Function to calculate overlap
    def calculate_overlap(seed_set, non_seed_set):
        return list(set(seed_set) & set(non_seed_set))

    # Parallelized function to calculate overlap for one row in df1 with all rows in df2
    def calculate_overlap_parallel(row1, df2):
        return [calculate_overlap(row1['SeedSet'], row2['NonSeedSet']) for _, row2 in df2.iterrows()]

    # Create a new DataFrame for overlaps
    # NOTE: iterate over the subsets !!!!!!!!!
    overlaps_df = pd.DataFrame(index=df14_subset.index, columns=df2.index)

    # Parallelize the overlap calculation using joblib with tqdm for progress tracking
    num_cores = -1  # Use all available cores

    results = Parallel(n_jobs=num_cores)(
        delayed(calculate_overlap_parallel)(row1, df2)
        for _, row1 in tqdm(df14_subset.iterrows(),
                            total=len(df14_subset)
                            )
        )

    # Fill in the overlaps DataFrame with calculated values
    for i, row in enumerate(results):
        overlaps_df.iloc[i] = row

    with open("overlaps_14.pckl","wb") as f:
        pickle.dump(overlaps_df, f)

    del results


# overlaps_df = pd.concat([overlaps_df_1, overlaps_df_2], ignore_index=True)
# print(overlaps_df)


elif sys.argv[1] == "format":

    """
    PART D:

    """
    import pickle
    import pandas as pd

    with open("overlaps_6.pckl","rb") as f:
        df = pickle.load(f)

    df = df.rename_axis(index={'PATRIC': 'PATRIC_beneficary'})
    df = df.rename_axis(columns={'PATRIC': 'PATRIC_donor'})

    stacked_df = df.stack()
    result_df = stacked_df.reset_index()
    result_df = result_df[['PATRIC_beneficary', 'PATRIC_donor']]
    # ATTENTION: Replace 1 with where we have left
    result_df['auto_increasing_number'] = range(307550001 + 1 , 246040001 + 1 + len(result_df))


    existing_csv_path = "seed_complements_map.tsv"
    result_df.to_csv(existing_csv_path, mode='a', sep="\t", header=False, index=False)


    # to get the 3-cols .tsv with the seeds
    output_file = "overlaps_14.tsv"
    df = df.rename_axis(index={'PATRIC': 'PATRIC_beneficary'})
    df = df.rename_axis(columns={'PATRIC': 'PATRIC_donor'})
    stacked_df = df.stack()
    result_df = stacked_df.reset_index()
    result_df.columns = ["PATRIC_beneficary", "PATRIC_donor", "seed_complement"]
    result_df.to_csv(output_file, sep='\t', index=False)



    f = open("seed_complements_map.tsv", "w")
    f.write("PATRIC_beneficary\tPATRIC_donor\tseedComplementId\n")
    all_ids = list(df.columns)
    counter = 0
    for i in all_ids:
        for j in all_ids:
            counter += 1
            f.write(i + "\t" + j + "\t" + str(counter) + "\n")

else:
    print(mhelp)




# awk 'BEGIN {OFS="\t"} {print $1, $2, ++count}' a.tsv > output.tsv
# awk -F"\t"  'NR==FNR { map[$2]=$1; next } { if ($1 in map) $1=map[$1]; if ($2 in map) $2=map[$2]; print $1, $2, $3 }' gtdb2patricIds.tsv overlaps_07.tsv > gtdb_overlaps_07.tsv