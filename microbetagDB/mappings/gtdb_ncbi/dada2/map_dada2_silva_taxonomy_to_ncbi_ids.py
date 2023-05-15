#!/usr/bin/env python
import pandas as pd

# ncbi_names_ids = pd.read_csv("/home/luna.kuleuven.be/u0156635/github_repos/microbetag/mappings/gtdb_ncbi/ncbi_dump/names.dmp", sep="\t|\t", header = None, engine="python")
# dada2_silva = open("silva_138_taxonomies_unique.tsv", "r")
# output = open("silva_138_taxonomies_ncbiIds.tsv", "w")
# non_match = open("non_matched_taxonomies.tsv", "w")

# for line in dada2_silva:
#     taxonomy = line.split(";")
#     check = False
#     if len(taxonomy) == 8:
#         taxon = " ".join([taxonomy[5], taxonomy[6]])
                
#     else: 
#         taxon = taxonomy[-2]

#     try:
#         corresponding_ncbi_id = ncbi_names_ids.loc[ ncbi_names_ids[2] == taxon ]
#     except KeyError:
#         print("No NCBI Id for species: ", taxon)


#     if corresponding_ncbi_id.empty:
#         non_match.write(line + "\n")
#         print(taxon)

#     else:
#         ncbi_id = corresponding_ncbi_id[0].values.tolist()[0]
        
#         # lines has a new line character in the end
#         output.write("\t".join([line[:-1], str(ncbi_id)]) + "\n")


# try again
remaining_taxonomies = open("non_matched_taxonomies.tsv", "r")
gtdb_metadata        = pd.read_csv("/home/luna.kuleuven.be/u0156635/github_repos/microbetag/mappings/gtdb_ncbi/gtdb_metadata/GTDB_QUALITY_REPRESENTATIVE_GENOMES_v202", sep="\t", header=None, low_memory=False)
species_left = open("species_left.tsv", "w")
tmp = open("tmp", "w")
for line in remaining_taxonomies: 


    taxonomy = line.split(";")

    if len(taxonomy) == 7:

        if "-" in taxonomy[5]: 

            taxon = " ".join([taxonomy[5].split("-")[-1], taxonomy[6]])

        else:

            taxon = " ".join([taxonomy[5], taxonomy[6]])

        species_left.write(taxon)
        tmp.write(line)

    else:

        print(taxonomy)

# sed -i 's/Candidatus //g' species_left.tsv

# while IFS= read -r line;  do grep  ";$line"  ../gtdb_metadata/GTDB_QUALITY_REPRESENTATIVE_GENOMES_v202 | awk -F"\t" -v mytaxon="$line" '{print mytaxon "\t" $15 "\t" $78}'  ; done < species_left.tsv
