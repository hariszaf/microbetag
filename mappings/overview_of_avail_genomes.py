#!/usr/bin/env python3

"""
This script will return a .json file (microbetag_ko_db.json) including 
the NCBI Taxonomy Ids that are currently available 
in the microbetag framework along with the the followin
information about each NCBI Id: 
- species or strain 
- resources from which a genome is available, currently: MGnify catalogues, KEGG GENOMES, GTDB 
"""
import os 
import re
import json
import sys

""" 
Keep track of the metadata related to the MGnify catalogues' MAGs retrieved
"""
cwd = os.getcwd()

mgnify_mag_metadata = {}
catalogues_path     = os.path.join(cwd, "mgnify_mappings/metadata")
catalogues          = os.listdir(catalogues_path)

for metadata_filemame in catalogues: 

    catalogue = metadata_filemame.split("_v")[0]
    version   = re.findall("_v.*_", metadata_filemame)

    catalogue_file = os.path.join(catalogues_path, metadata_filemame)

    for line in open(catalogue_file, "r"): 

        line = line.split("\t")
        mag  = line[0]
        link = line[19]

        mgnify_mag_metadata[mag]              =  {}
        mgnify_mag_metadata[mag]['catalogue'] = catalogue
        mgnify_mag_metadata[mag]['version']   = version
        mgnify_mag_metadata[mag]['link']      = link

"""
Keep track of the taxonomic level of the NCBI Taxonomy ids found
"""

taxdump           = open("taxdump/nodes.dmp", "r")
ncbi_id_tax_level = {}

for line in taxdump:

    line      = line.split("\t|\t")
    tax_id    = line[0]
    tax_level = line[2]

    ncbi_id_tax_level[tax_id] = tax_level



"""
Parse the genomes retrieved 
"""
ncbi_ids_tax_level_resources = {}
list_of_resources            = ["kegg_genomes", "mgnify_catalogues"]    # once annotations completed, add "gtdb_genomes/"
mgnify_url_base              = "http://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/"
for resource in list_of_resources:


    for ncbi_id in os.listdir(resource): 

        try: 

            if ncbi_id not in ncbi_ids_tax_level_resources: 
                ncbi_ids_tax_level_resources[ncbi_id]                                  = {}
                ncbi_ids_tax_level_resources[ncbi_id]["tax-level"]                     = ncbi_id_tax_level[ncbi_id]
                ncbi_ids_tax_level_resources[ncbi_id]["resources"]                     = [resource]
                ncbi_ids_tax_level_resources[ncbi_id]["links"]                         = {}
                ncbi_ids_tax_level_resources[ncbi_id]["links"][resource]               = {}

            else:

                if resource not in ncbi_ids_tax_level_resources[ncbi_id]["resources"]:

                    ncbi_ids_tax_level_resources[ncbi_id]["resources"].append(resource)                    
                    ncbi_ids_tax_level_resources[ncbi_id]["links"][resource]               = {}

        except: 

            print("step 1: NCBI Taxonomy id abbandoned ", ncbi_id)


        for filename in os.listdir(os.path.join(resource, ncbi_id)):

            if "kegg" in resource:
                kegg_ids   = []

                if filename.endswith('_kos_related_to_mos.json'):
                    kegg_ids.append(filename.split("_kos_related_to_mos.json")[0])

                for kegg_id in kegg_ids: 
                    kegg_file  = kegg_id + "_kos_related_to_mos.json"
                    inner_link = os.path.join(cwd, "ref-dbs/kegg_genomes/", ncbi_id, kegg_file)
                    outer_link = "https://rest.kegg.jp/link/ko/" + kegg_id

                    ncbi_ids_tax_level_resources[ncbi_id]["links"][resource][kegg_id] = {}
                    ncbi_ids_tax_level_resources[ncbi_id]["links"][resource][kegg_id]['inner-link'] = kegg_file
                    ncbi_ids_tax_level_resources[ncbi_id]["links"][resource][kegg_id]['outer-link'] = outer_link



            if "mgnify" in resource: 

                mags_ids   = []
                if filename.endswith('_kos_related_to_mos.json'): 
                    mags_ids.append(filename.split('_kos_related_to_mos.json')[0])

                for mag_id in mags_ids: 

                    mag_file   = mag_id + "_kos_related_to_mos.json"

                    inner_file = os.path.join(cwd, "ref-dbs/mgnify_catalogues/", ncbi_id, mag_file)

                    metadata_link     = mgnify_mag_metadata[mag_id]['link']
                    catalogue         = metadata_link.split("/")[7]
                    catalogue_version = metadata_link.split("/")[8]
                    outer_link        = os.path.join(mgnify_url_base, catalogue, catalogue_version, "species_catalogue", mag_id[:-2], mag_id, "genome")

                    try:
                        ncbi_ids_tax_level_resources[ncbi_id]["links"][resource][mag_id] = {}
                        ncbi_ids_tax_level_resources[ncbi_id]["links"][resource][mag_id]['inner-link'] = mag_file
                        ncbi_ids_tax_level_resources[ncbi_id]["links"][resource][mag_id]['outer-link'] = outer_link
                    except:
                        print("step 2: NCBI Taxonomy id abbandoned ", ncbi_id)

with open("microbetag_ko_db.json", "w") as json_file:
    json.dump(ncbi_ids_tax_level_resources, json_file, indent = 4)



# """
# if the NCBI Taxonomy id corresponds to a species, then it is the UNION 
# of all the genomes found that need to be taken into account, meaning the KO terms found. 
# if the Id stands for a strain, then it is the INTERSECTION of the 
# KO terms that will be considered. 
# """

# for ncbi_id in ncbi_ids_tax_level_resources.items():

#     if ncbi_ids_tax_level_resources[ncbi_id]['tax-level'] == "species": 



#     else: 



