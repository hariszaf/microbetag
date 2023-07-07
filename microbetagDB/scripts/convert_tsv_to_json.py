
import glob 
import os, sys
import json
import shutil

all_genomes_dir = "/home/luna.kuleuven.be/u0156635/github_repos/microbetag/ref-dbs/ALL_GENOME_MODULES"

"""
convert .tsv files from gtdb_genomes folder to .json and move them under ALL_GENOME_MODULES
move the already .json files under kegg_genomes to ALL_GENOME_MODULES
likewise, for the mgnify_catalogues case
"""

# for gdir in os.listdir("gtdb_genomes/"): 

#     gpath = os.path.join(all_genomes_dir, gdir) 

#     if not os.path.exists(gpath):
#         os.mkdir(gpath)

#     for gfile in glob.glob("gtdb_genomes/" + gdir + "/" + "*_related_to_mos.tsv"):

#         genome_name = gfile.split("/")[-1].split(".")[0]
#         genome_name = genome_name + "." + gfile.split("/")[-1].split(".")[1].split("_")[0]

#         f = open(gfile, "r")
#         g = f.readlines()
#         gdic = {}

#         for pair in g:
#             g = pair.split("\t")
#             md = g[0]
#             term = g[1][:-1]

#             if md not in gdic:
#                 gdic[md] = [term]

#             else:
#                 gdic[md].append(term)
    
#         true_gdic = {}
#         for gmodule, gterms in gdic.items():

#             counter = 0
#             pif     = {}
#             for term in gterms:
#                 pif[counter] = term 
#                 counter += 1


#             true_gdic[gmodule] = pif


#         ooout = "".join([genome_name, "_kos_related_to_mos.json"])

#         with open(os.path.join(gpath, ooout), "w") as outfile:
#             json.dump(true_gdic, outfile)

# can be used for the kegg_genomes/ case as well
for gdir in os.listdir("mgnify_catalogues"): 

    gpath = os.path.join(all_genomes_dir, gdir) 

    if not os.path.exists(gpath):
        os.mkdir(gpath)

    for gfile in glob.glob("mgnify_catalogues/" + gdir + "/" + "*_related_to_mos.json"):

        copy_gfile = os.path.join(gpath, gfile.split("/")[-1])
        shutil.copyfile(gfile, copy_gfile)
