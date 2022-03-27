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

cwd             = os.getcwd()
microbetag_path = os.path.dirname(cwd)
db_path         = os.path.join(microbetag_path, "ref-dbs/")


""" 
Keep track of the metadata related to the MGnify catalogues' MAGs retrieved
"""

mgnify_mag_metadata = {}
catalogues_path     = os.path.join(cwd, "mgnify_mappings/metadata")
catalogues          = os.listdir(catalogues_path)

for metadata_filemame in catalogues: 

    catalogue = metadata_filemame.split("_v")[0]
    version   = re.findall("_v.*_", metadata_filemame)

    catalogue_file = os.path.join(catalogues_path, metadata_filemame)
    with open(catalogue_file, "r") as f:
        next(f)
        for line in f: 

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
taxdump_file      = os.path.join(db_path, "taxdump/nodes.dmp")
taxdump           = open(taxdump_file, "r")
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

    for ncbi_id in os.listdir(os.path.join(db_path, resource)): 

        try: 

            if ncbi_id not in ncbi_ids_tax_level_resources: 
                ncbi_ids_tax_level_resources[ncbi_id]                    = {}
                ncbi_ids_tax_level_resources[ncbi_id]["tax-level"]       = ncbi_id_tax_level[ncbi_id]
                ncbi_ids_tax_level_resources[ncbi_id]["resources"]       = [resource]
                ncbi_ids_tax_level_resources[ncbi_id]["links"]           = {}
                ncbi_ids_tax_level_resources[ncbi_id]["links"][resource] = {}

            else:

                if resource not in ncbi_ids_tax_level_resources[ncbi_id]["resources"]:

                    ncbi_ids_tax_level_resources[ncbi_id]["resources"].append(resource)                    
                    ncbi_ids_tax_level_resources[ncbi_id]["links"][resource] = {}

        except: 
            print("step 1: NCBI Taxonomy id abbandoned ", ncbi_id)


        for filename in os.listdir(os.path.join(db_path, resource, ncbi_id)):

            if "kegg" in resource:
                kegg_ids   = []

                if filename.endswith('_kos_related_to_mos.json'):
                    kegg_ids.append(filename.split("_kos_related_to_mos.json")[0])

                for kegg_id in kegg_ids: 
                    kegg_file  = kegg_id + "_kos_related_to_mos.json"
                    inner_link = os.path.join("ref-dbs/kegg_genomes/", ncbi_id, kegg_file)
                    outer_link = "https://rest.kegg.jp/link/ko/" + kegg_id


                    ncbi_ids_tax_level_resources[ncbi_id]["links"][resource][kegg_id] = {}
                    ncbi_ids_tax_level_resources[ncbi_id]["links"][resource][kegg_id]['inner-link'] = inner_link
                    ncbi_ids_tax_level_resources[ncbi_id]["links"][resource][kegg_id]['outer-link'] = outer_link

            if "mgnify" in resource: 

                mags_ids   = []
                if filename.endswith('_kos_related_to_mos.json'): 
                    mags_ids.append(filename.split('_kos_related_to_mos.json')[0])

                for mag_id in mags_ids: 

                    mag_file   = mag_id + "_kos_related_to_mos.json"
                    inner_file = os.path.join("ref-dbs/mgnify_catalogues/", ncbi_id, mag_file)

                    metadata_link     = mgnify_mag_metadata[mag_id]['link']
                    catalogue         = metadata_link.split("/")[7]
                    catalogue_version = metadata_link.split("/")[8]
                    outer_link        = os.path.join(mgnify_url_base, catalogue, catalogue_version, "species_catalogue", mag_id[:-2], mag_id, "genome")

                    try:
                        ncbi_ids_tax_level_resources[ncbi_id]["links"][resource][mag_id] = {}
                        ncbi_ids_tax_level_resources[ncbi_id]["links"][resource][mag_id]['inner-link'] = inner_file
                        ncbi_ids_tax_level_resources[ncbi_id]["links"][resource][mag_id]['outer-link'] = outer_link
                    except:
                        print("step 2: NCBI Taxonomy id abbandoned ", ncbi_id)

with open("microbetag_ko_db.json", "w") as json_file:
    json.dump(ncbi_ids_tax_level_resources, json_file, indent = 4)

"""
If the NCBI Taxonomy id corresponds to a species, then it is the UNION 
of all the genomes found that will be taken into account (by default), meaning all the KO terms found. 
If the Id stands for a strain, then it is only the INTERSECTION of the KO terms that will be considered. 

Here we make a file for each NCBI Taxonomy id retrieved, with all the KO terms 
for all its corresponding genomes (under the /ref-dbs/microbetag_genomes/union/ path)
and another one with only those terms that are present in all its corresponding genomes 
(under the /ref-dbs/microbetag_genomes/union/intersection path).
"""

for ncbi_id in ncbi_ids_tax_level_resources:


    genomes_avail = set()

    for category in ncbi_ids_tax_level_resources[ncbi_id]['links']:

        for genome in ncbi_ids_tax_level_resources[ncbi_id]['links'][category]: 

            try: 
                genomes_avail.add(ncbi_ids_tax_level_resources[ncbi_id]['links'][category][genome]['inner-link'])
            except:
                print(ncbi_ids_tax_level_resources[ncbi_id])

    # If a KO appears in the kos_found list as many times as the len(genomes_aval) then it belongs to the intersection
    # The unique elements of the kos_found list, make up the union case
    # Attention though! A KO term might appears in a single 'genome' multiple times, thus we need to count it only once per 'genome'
    kos_found        = []
    kos_md           = set()

    for genome in genomes_avail:

        kos_per_md    = json.load(open(os.path.join(microbetag_path, genome), "r"))
        current_kos   = set()

        for md, kos in kos_per_md.items(): 
            for ko in kos.values():
                current_kos.add(ko)
                kos_md.add((ko, md))

        kos_found.extend(current_kos)

    union_kos        = set()
    intersection_kos = set()

    for ko in kos_found: 
        if kos_found.count(ko) == len(genomes_avail): 
            intersection_kos.add(ko)
            union_kos.add(ko)
        else:
            union_kos.add(ko)

    union_kos_per_md        = {}
    intersection_kos_per_md = {}

    """ 
    Make a file for the union case.
    """
    for ko in union_kos: 

        for pair in kos_md: 

            if ko in pair: 

                md = pair[1]

                if md not in union_kos_per_md: 
                    union_kos_per_md[md]      = {}
                    union_kos_per_md[md]['0'] = ko

                else: 
                    union_kos_per_md[md][str(len(union_kos_per_md[md]))] = ko

    union_filename  = ncbi_id + "_" + "union_kos_related_to_mos.json"
    union_ncbi_path = os.path.join(db_path, "microbetag_genomes/union", ncbi_id)
    union_file      = os.path.join(db_path, "microbetag_genomes/union", ncbi_id, union_filename)

    if os.path.isdir(union_ncbi_path) == False:
        os.mkdir(union_ncbi_path)

    with open(union_file, "w") as g: 
        json.dump(union_kos_per_md, g, indent=4)


    """
    Make a file for the intersection case.  
    """
    for ko in intersection_kos: 
        for pair in kos_md: 
            if ko in pair: 
                md = pair[1]
                if md not in intersection_kos_per_md: 
                    intersection_kos_per_md[md]      = {}
                    intersection_kos_per_md[md]['0'] = ko

                else: 
                    intersection_kos_per_md[md][str(len(intersection_kos_per_md[md]))] = ko

    intersection_filename  = ncbi_id + "_" + "intersection_kos_related_to_mos.json"
    intersection_ncbi_path = os.path.join(db_path, "microbetag_genomes/intersection", ncbi_id)
    intersection_file      = os.path.join(db_path, "microbetag_genomes/intersection", ncbi_id, intersection_filename)

    if os.path.isdir(intersection_ncbi_path) == False:
        os.mkdir(intersection_ncbi_path)

    with open(intersection_file, "w") as g: 
        json.dump(intersection_kos_per_md, g, indent=4)


