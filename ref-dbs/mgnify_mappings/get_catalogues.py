#!/usr/bin/env python
"""
In case you want to use the MGnify JSON:API, have a loot at:
https://github.com/EBI-Metagenomics/examples/blob/master/mgnify/src/notebooks/answers/ANSWER_examples.ipynb

The jsonapi_client (https://github.com/qvantel/jsonapi-client) is a library
that enables the use of a JSON
"""

import os, json, time, requests
import random
from jsonapi_client import Session # pip install jsonapi-client
import logging
from subprocess import Popen, PIPE

# Link biosample or ncbi genome id to ncbi taxonomy id
def efetch_using_ncbi_genome_id(genome_id, db):
   """
   This function is no longer needed. 
   It is part of the get_ncbi_tax_id() function that has been replaced
   """

   try:

      logging.info("About to get the NCBI Taxonomy Id of the corresponding species")

      p1 = Popen(['esearch', '-db', db, '-query', genome_id], stdout=PIPE)
      p2 = Popen(['elink', '-target', 'taxonomy'],               stdin=p1.stdout, stdout=PIPE)
      p3 = Popen(['esummary'], stdin=p2.stdout, stdout=PIPE)
      p4 = Popen(['xtract', '-pattern', 'DocumentSummary', '-element', 'TaxId'],  stdin=p3.stdout, stdout=PIPE)
      stdout, _ = p4.communicate()
      ncbi_tax_id = stdout.decode("utf-8")[:-1]
      time.sleep(random.uniform(0.05, 0.2))
      print(ncbi_tax_id, type(ncbi_tax_id))
      return ncbi_tax_id

   except:
      logging.exception("Something wrong with the efetch function!")

# Get MAGs metadata
def get_quality_mgnify_mags(completeness_score, contamination_score):

   API_BASE = "https://www.ebi.ac.uk/metagenomics/api/latest/"

   count_non_quality_mags = 0
   quality_mgnify_genomes = {}

   with Session(API_BASE) as s:
      all_genomes_overall     = s.get('genomes')
      total_number_of_genomes = all_genomes_overall.meta.pagination['count']
      total_number_of_pages   = all_genomes_overall.meta.pagination['pages']

      for page in range(1, total_number_of_pages + 1):
         print("\n\n\n ~~~~~~ Page: ", str(page), "out of ", str(total_number_of_pages))

         with Session(API_BASE) as t:
            page_url = 'genomes?page=' + str(page)
            all_page_genomes_overall     = t.get(page_url)
            time.sleep(0.4)

            for genome in all_page_genomes_overall.resources: 
               if genome.completeness > 90.0 and genome.contamination < 5.0:
                  quality_mgnify_genomes[genome.accession] = {} # MGnigy genome accesion number as id
                  quality_mgnify_genomes[genome.accession]['ncbi_genome_accession_number'] = genome.ncbi_genome_accession
                  print(genome.ncbi_genome_accession)
                  quality_mgnify_genomes[genome.accession]['ena_sample_accession']         = genome.ena_sample_accession
                  quality_mgnify_genomes[genome.accession]['catalogue']                    = genome.catalogue.as_resource_identifier_dict()['id'] 
                  quality_mgnify_genomes[genome.accession]['lineage']                      = genome.taxon_lineage
               else:
                  count_non_quality_mags += 1
   print("Non quality mags: ", count_non_quality_mags, " out of: ", total_number_of_genomes)
   return quality_mgnify_genomes

# Get NCBI Taxonomy Id
def get_ncbi_tax_id(metadata_dict):

   """
   This function is no longer needed. 
   An alternative way will be exploited. See function map_gtdb_lineage_to_ncbi_tax_id()
   """

   for genome, metadata in metadata_dict.items():

      gca = str(metadata_dict[genome]['ncbi_genome_accession_number'])

      ena = str(metadata_dict[genome]['ena_sample_accession'])

      if gca != "null" and gca != "None":
         
         try:
            ncbi_tax_id = efetch_using_ncbi_genome_id(gca, "assembly")

            metadata_dict[genome]['ncbi_tax_id'] = ncbi_tax_id

         except:
            logging.exception("Something wrong with the efetch function!")


      elif ena != "null" and ena != "None":

         try:
            ncbi_tax_id = efetch_using_ncbi_genome_id(ena, "biosample")

            metadata_dict[genome]['ncbi_tax_id'] = ncbi_tax_id

         except:
            logging.exception("Something wrong with the efetch function!")

   return metadata_dict

# Get NCBI Taxonomy Id for the MAGs with GTDB lineage at the species level
def map_gtdb_lineage_to_ncbi_tax_id(metadata_dict):
   """
   This is valid only if you are interested for GTDB lineages at the species level
   In case that you would like a higher level taxonomy group, you cannot use this function
   """

   gtdb_lineage_ncbi_id = {}
   gtdb_metadata_file   = open("gtdb_lineage2ncbi_tax_id.tsv", "r")
   next(gtdb_metadata_file)

   for line in gtdb_metadata_file:
      gtdb_lineage_ncbi_id[line.split("\t")[0]] = line.split("\t")[1][:-1]

   for genome, metadata in metadata_dict.items():

      lineage = str(metadata_dict[genome]['lineage'])

      if lineage in gtdb_lineage_ncbi_id.keys():

         metadata_dict[genome]['ncbi_tax_id'] = gtdb_lineage_ncbi_id[lineage]

      else:

         metadata_dict[genome]['ncbi_tax_id'] = ""

   return(metadata_dict)

# Get the KOs assigned to the MAG
def get_KOs(metadata_dict): 

   for mag, metadata in metadata_dict.items():
      # print(mag)
      API_BASE = "http://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/"
      
      catalogue = "-".join(metadata_dict[mag]['catalogue'].split("-")[:-2])

      version   = ".".join(metadata_dict[mag]['catalogue'].split("-")[-2:])  
      
      url = API_BASE + catalogue + "/" + version + "/species_catalogue/" + mag[:-2] + "/" + mag + "/" + "genome/" + mag + "_eggNOG.tsv"
      
      r = requests.get(url)
      print(r.text)

      sys.exit(0)




## STEP A
# genomes_metadata = get_quality_mgnify_mags(90, 5)
# print("\n\n", genomes_metadata)
# initial_mag_dic = open("genomes_metadata.json", "w")  
# json.dump(genomes_metadata, initial_mag_dic, indent = 6)
# initial_mag_dic.close()


## STEP B
init_file       = open("genomes_metadata.json", "r")
init_dic        = json.load(init_file)

mags_with_ncbi_ids = map_gtdb_lineage_to_ncbi_tax_id(init_dic)


get_KOs(mags_with_ncbi_ids)


















"""
if NCBI Genome Accession available:
esearch -db assembly -query "GCA_015260835" | elink -target taxonomy | esummary | xtract -pattern DocumentSummary -element TaxId 

is ENA sample accession available: 
esearch -db biosample -query "ERS7766781" | elink -target taxonomy | esummary | xtract -pattern DocumentSummary -element TaxId 
"""



#  NOT TO USE !!!! 


# from IPython.display import Markdown
# import pandas as pd
# def get_variable_from_link_or_input(variable, name = 'accession', default = None):
#     """
#     Get a variable value, either from an ENV VAR that would have been set by the shiny_proxy_jlab_query_parms extension, or through direct user input.
#     """
#     var = os.getenv(variable)
#     if var:
#         display(Markdown(f'<span style="background-color: #0a5032; color: #fff; padding: 8px;">Using {name} <emph>{var}</emph> from the link you followed.</span>'))
#     else:
#         var = input(f'Type {"an" if name[0].lower() in "aeiou" else "a"} {name} [default: {str(default)}]')
#     var = var or default
#     print(f'Using "{var}" as {name}')
#     return var


# # api_endpoint = get_variable_from_link_or_input('MGYS00000596', 'API Endpoint', 'super-studies')
# accession = get_variable_from_link_or_input('MGYS', 'Study Accession', 'MGYS00005292')


# # # with Session("https://www.ebi.ac.uk/metagenomics/api/v1") as mgnify:
# # #     resources = map(lambda r: r.json, mgnify.iterate(api_endpoint))
# # #     resources = pd.json_normalize(resources)
# # #     resources.to_csv(f"{api_endpoint}.csv")
# # # resources

# df = DataFrame(columns=('category', 'description', 'annotation counts'))
# df.index.name = 'GO term'

# with Session(API_BASE) as s:
#     run = s.get('runs', 'ERR598955').resource
#     print("Run: ", run, "\t type of run: ", type(run))
#     for a in run.analyses:
#         for ann in a.go_slim:
#             df.loc[ann.accession] = [
#                 ann.lineage, ann.description, ann.count
#             ]
# df