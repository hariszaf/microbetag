#!/usr/bin/env python3

__author__  = "Haris Zafeiropoulos"
__email__   = "haris.zafeiropoulos@kuleuven.be"

"""
Aim: 

This script makes use of the ncbi-taxonomist tool to map the GTDB taxa to their corresponding NCBI Taxonomy Ids if such exist. 

Usage:

./ncbi_taxonomist_for_genus.py <one_column_file_with_gtdb_taxa_names> <output_filename>

Requires: 

* the ncbi-taxonomist library that handles and manages phylogenetic data from NCBI's Entrez taxonomy database (https://gitlab.com/janpb/ncbi-taxonomist)
        pip install ncbi-taxonomist

* the jq package; a lightweight and flexible command-line JSON processor. To get it:
        sudo apt-get install jq
"""

import subprocess, time
import sys
import time 

taxonomies_file = open(sys.argv[1], "r")
input_names     = taxonomies_file.readlines()
output_file     = open(sys.argv[2], "w")


for line in input_names:
 
    check = False
    if "__" in line:
        gtdb_taxonomy  = line.split("__")[1]
    else: 
        gtdb_taxonomy = line

    print(gtdb_taxonomy)

    ps     = subprocess.Popen(('ncbi-taxonomist', 'resolve', '-n', gtdb_taxonomy), stdout=subprocess.PIPE)
    output = subprocess.check_output(("jq", ".taxon.taxid"), stdin=ps.stdout, stderr=subprocess.STDOUT)
    ps.wait()

    time.sleep(0.12)

    output = str(output)
    output = output[2:]
    output = output[:-1]

    if output != "":

        output = output[:-2]

        if " " not in output:

            ncbi_id = output 
            check = True
            output_file.write(line[:-1] + "\t" + ncbi_id + "\n")

    if not check:
        output_file.write(line[:-1] + "\t" + "nan" + "\n")
