#!/usr/bin/env python
"""
This script will download all the corresponding .faa files 
to the high quality (completeness > 95% and contamination < 5%) representative genomes of GTDB. 

Then, the .faa files will be annotated for KEGG terms using HMMER (http://hmmer.org/)
"""

import os
import time
import sys
import requests

metadata_file = open("GTDB_QUALITY_REPRESENTATIVE_GENOMES", "r")

API_BASE = "https://ftp.ncbi.nlm.nih.gov/genomes/all/"

counter = 0
for line in metadata_file:

   counter +=1 

   if counter % 10 == 0:
      print("Genomes parsed: ", str(counter), "out of 26,778")

   if counter <= 20557:
       continue

   ncbi_genome_accession = line.split("\t")[54]
   ncbi_assembly         = line.split("\t")[46]

   url_letters         = ncbi_genome_accession.split("_")[0]
   numeric_part        = ncbi_genome_accession.split("_")[1].split(".")[0]

   split_to_num_parts  = [numeric_part[i:i+3] for i in range(0, len(numeric_part), 3)]
   url_numbers         = '/'.join(split_to_num_parts)

   url = API_BASE + url_letters + "/" + url_numbers + "/" + ncbi_genome_accession + "_" + ncbi_assembly + "/" + ncbi_genome_accession + "_" + ncbi_assembly + "_protein.faa.gz"

   ref_dbs_dir = os.path.dirname(os.getcwd())
   assembly_saving_name = ncbi_assembly.replace("/", "_")
   faa_file    = ref_dbs_dir + "/gtdb_genomes/faa_files/" + ncbi_genome_accession + "_" + assembly_saving_name + ".faa.gz"
   r = requests.get(url, allow_redirects = True)
   open(faa_file, 'wb').write(r.content)

   time.sleep(0.5)



