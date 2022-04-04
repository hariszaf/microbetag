#!/usr/bin/env python3

import subprocess, time
import sys
import time 

# genus_names   = open("genus.uniq", "r")
family_names  = open("family.uniq.no_genus", "r")
level_ncbi_id = open("family_names_to_ncbiId.tsv", "w")

input_names = family_names.readlines()

for line in input_names:
 
    level  = line.split("__")[1]

    ps     = subprocess.Popen(('/home/haris/anaconda3/bin/ncbi-taxonomist', 'resolve', '-n', level), stdout=subprocess.PIPE)
    output = subprocess.check_output(("jq", ".taxon.taxid"), stdin=ps.stdout, stderr=subprocess.STDOUT)
    ps.wait()

    time.sleep(0.2)

    output = str(output)
    output = output[2:]
    output = output[:-1]

    print(output)

    if output != "":

        output = output[:-2]

        if " " not in output:

            ncbi_id = output 
            level_ncbi_id.write(line[:-1] + "\t" + ncbi_id + "\n")

        print(level + "\t" + ncbi_id)

    print("\n\n~~~\n\n")