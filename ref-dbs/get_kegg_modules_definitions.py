#!/usr/bin/env python3

import time, sys
import requests

all_modules_ids     = []
# modules_definitions = {}

for line in open("kegg_terms_per_module.tsv", "r"): 
   md_id = line.split("\t")[0]
   if md_id not in all_modules_ids: 
      all_modules_ids.append(md_id)


print("# of modules: ", len(all_modules_ids))

for md in all_modules_ids:

   print(md)
   
   url       = "http://rest.kegg.jp/get/" + md
   response   = requests.get(url)
   # print(response.content)

   with open("response.txt", "wb") as f:
      f.write(response.content)
      
   for line in open("response.txt", "r"):
      if "DEFINITION" in line:
         definition = line.split("  ")[1:][0][:-1]
         definition = ';'.join(definition.split(" "))
         # modules_definitions[md] = definition

         with open("module_definitions.tsv", "a") as f:
               f.write(md + "\t" + definition + "\n")
   
   time.sleep(0.5)
         


# Definition	The definition of the module as a logcail expression of K numbers. 
# Comma separated K numbers indicate alternatives. 
# Plus signs are used to represent a complex and a minus sign denotes a non-essential component in the complex.



