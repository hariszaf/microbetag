#!/usr/bin/env python3

import glob
import json
import sys

def parse_definitions_file():

   definitions_file = open("../ref-dbs/module_definitions.tsv", "r")
   module_definitions_steps = {}
   

   for line in definitions_file:

      md         = line.split("\t")[0]
      definition = line.split("\t")[1][:-1]
      list_of_lists_of_single_steps = []

      if "(" in definition or "," in definition:

         print(md, definition)

         quasi_steps = definition.split(";")
         print("quasi steps:")
         print(quasi_steps)


         for quasi_step in quasi_steps:

            num_of_open   = quasi_step.count('(')
            num_of_closed = quasi_step.count(')')

            if num_of_open == num_of_closed:
               """
               every quasi step is an actual step
               this is so, cause if we have an inner step in a choice
               we will have an inequality in all cases
               """

               if num_of_open == 0:
                  list_of_lists_of_single_steps.append(quasi_step)

               else:

                  print("EDO KOITA MPINE")
                  print(quasi_step)
                  quasi_step = quasi_step.replace("(","")
                  quasi_step = quasi_step.replace(")","")
                  quasi_step = quasi_step.split(",")
                  
                  print(quasi_step)





         print("\n~~~~\n")







      else:
         """
         no alternatives per step 
         REMEMBER! A step can be a complex, so more than 1 KOs may be required for a single sptep, 
         meaning that a single step (a list) may include multiple lists 
         """
         steps = definition.split(";")
         for step in steps:
         
            # This denotes that there are reactions for whom there are no corresponding KO terms 
            # Check this if after discussion with KF
            if "--" in step:
               continue

            if "+" not in step and "-" not in step:
               list_of_lists_of_single_steps.append([step])

            elif "+" in step and "-" not in step:
               step = step.split("+")
               list_of_semis = [[semi] for semi in step]
               list_of_lists_of_single_steps.append(list_of_semis)

            elif "-" in step and "+" not in step: 

               step  = step.split("-")

               # A term you can ignore
               if len(step) == 2 and step[0] == "":
                  continue

               else:
                  list_of_lists_of_single_steps.append(step[0])

            # Check this if after discussion with KF
            elif "-" in step and "+" in step:
               
               semis = []
               step  = step.split("-")
               semi_counter = 0
               
               if len(step) == 2 and step[0] == "":
                  print("Ignore this optional KO term, even if there's a + ??")
                  continue

               elif len(step) == 2 and step[0] != "":
                  list_of_lists_of_single_steps.append(step[0])

               else: 
                  # step has more than 2 parts
                  for semi in step:
                     semi_counter += 1

                     if "+" in semi:

                        semi = semi.split("+")

                        if semi_counter == 1:
                           list_of_semis = [part for part in semi]
                           semis.append(list_of_semis)

                        else:
                           list_of_semis = [part for part in semi[1:]]
                           semis.append(list_of_semis)


                  filtered_step = []
                  for c in range(len(semis)):                     
                     for v in range(len(semis[c])):
                        filtered_step.append(semis[c][v])

                  list_of_lists_of_single_steps.append(filtered_step)


         y                                          = 0
         num_of_steps                               = len([y+1 for x in list_of_lists_of_single_steps if "--" not in x ])
         list_of_lists_of_single_steps              = tuple(list_of_lists_of_single_steps)
         module_definitions_steps[md]               = {}
         module_definitions_steps[md]['steps']      = list_of_lists_of_single_steps
         module_definitions_steps[md]['# of steps'] = num_of_steps




   # print(module_definitions_steps)
   with open("test.json", "w") as f:
      json.dump(module_definitions_steps, f)


def map_genome_to_modules(ncbi_taxonomy_id):

   ncbi_taxonomy_id = str(ncbi_taxonomy_id)

   module_definitions = open("../ref-dbs/module_definitions.tsv")
   kegg_genome_file   = glob.glob("../ref-dbs/kegg_genomes/" + ncbi_taxonomy_id + "/*.json") 
   kegg_genome_file   = open(kegg_genome_file[0])
   kegg_genome        = json.load(kegg_genome_file)



parse_definitions_file()
# map_genome_to_modules(36033)


