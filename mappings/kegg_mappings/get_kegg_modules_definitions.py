#!/usr/bin/env python3

import time, sys
import requests


def write_definition(line, md, iter):

   if iter > 1: 
      definition = line.strip()
      print("HERE WE GO:  ", definition)
      with open("module_definitions.tsv", "a") as f:
            f.write(md + "_" + str(iter) + "\t" + definition + "\n")

   else:
      definition = line.split("  ")[1:][0][:-1]
      definition = ';'.join(definition.split(" "))
      with open("module_definitions.tsv", "a") as f:
            f.write(md + "\t" + definition + "\n")

   return True



all_modules_ids     = []

for line in open("kegg_terms_per_module.tsv", "r"): 
   md_id = line.split("\t")[0]
   if md_id not in all_modules_ids: 
      all_modules_ids.append(md_id)

print("# of modules: ", len(all_modules_ids))

for md in all_modules_ids:
   
   print(md)

   url       = "http://rest.kegg.jp/get/" + md
   response   = requests.get(url)

   with open("response.txt", "wb") as f:
      f.write(response.content)

   switch      = False
   num_of_defs = 0 

   for line in open("response.txt", "r"):

      if "DEFINITION" in line:
         switch      = True
         num_of_defs = 1
         write_definition(line, md, num_of_defs)
         continue

      if "ORTHOLOGY" in line: 
         switch = False
      
      if switch: 
         num_of_defs += 1
         write_definition(line, md, num_of_defs)
   
   time.sleep(0.5)
         



# Definition	The definition of the module as a logcail expression of K numbers. 
# Comma separated K numbers indicate alternatives. 
# Plus signs are used to represent a complex and a minus sign denotes a non-essential component in the complex.



