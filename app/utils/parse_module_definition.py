#!/usr/bin/env python3

import glob
import json
import sys
from itertools import product


def parse_definitions_file():

   definitions_file = open("../ref-dbs/module_definitions.tsv", "r")
   module_definitions_steps = {}
   

   for line in definitions_file:

      md           = line.split("\t")[0]
      definition   = line.split("\t")[1][:-1]
      parsed_steps = []



      if "(" not in definition and "," not in definition:

         """
         no alternatives per step 
         REMEMBER! A step can be a complex, so more than 1 KOs may be required for a single sptep, 
         meaning that a single step (a list) may include multiple lists 
         """
         steps = definition.split(";")
         parsed_steps = parse_valid_steps_of_a_module(steps)


      else:

         intermediate_steps = []

         quasi_steps = definition.split(";")

         # print("\n\n\n")
         # print(md,definition)
         # print("------------------")

         count_front_parenth   = 0
         count_reverse_parenth = 0
         
         pseudo_step = []
         
         for quasi_step in quasi_steps:

            in_parenthesis = False

            if "(" not in quasi_step and ")" not in quasi_step and in_parenthesis == False:
               
                  intermediate_steps.append(parse_valid_steps_of_a_module([quasi_step]))

            else:

               if "(" in quasi_step or ")" in quasi_step:

                  num_of_open            = quasi_step.count('(')
                  count_front_parenth   += num_of_open
                  num_of_closed          = quasi_step.count(')')
                  count_reverse_parenth += num_of_closed
                  in_parenthesis         = True

               if count_front_parenth > count_reverse_parenth:
                  pseudo_step.append(quasi_step)
               
               else:
                  pseudo_step.append(quasi_step)
                  # intermediate_steps.append(parse_valid_steps_of_a_module(pseudo_step))

                  complete_step = ';'.join(pseudo_step) 


                  # print("definition: ", definition)
                  # print("quasi step: ", quasi_step)

                  # print("COMPLEX STEP: ", complete_step)

                  break_down_complex_step(complete_step, definition, md)



                  pseudo_step           = []
                  count_front_parenth   = 0
                  count_reverse_parenth = 0
                  in_parenthesis        = False


      y                                          = 0
      num_of_steps                               = len([y+1 for x in parsed_steps if "--" not in x ])
      parsed_steps              = tuple(parsed_steps)
      module_definitions_steps[md]               = {}
      module_definitions_steps[md]['steps']      = parsed_steps
      module_definitions_steps[md]['# of steps'] = num_of_steps


   # print(module_definitions_steps)
   with open("test.json", "w") as f:
      json.dump(module_definitions_steps, f)



def parse_valid_steps_of_a_module(steps):

   list_of_lists_of_single_steps = []

   for step in steps:

      # This denotes that there are reactions for whom there are no corresponding KO terms 
      # Check this if after discussion with KF
      if "--" in step:
         continue

      if "+" not in step and "-" not in step:
         list_of_lists_of_single_steps.append([step])

      elif "+" in step and "-" not in step:
         step = step.split("+")
         list_of_semis = [semi for semi in step]
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
         
         # print(step)
         semis = []
         step  = step.split("-")
         semi_counter = 0
         
         if len(step) == 2 and step[0] != "":
            if "+" in step[0]:
               components = step[0].split("+")
               list_of_lists_of_single_steps.append(components)


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
            # print(filtered_step)
            list_of_lists_of_single_steps.append(filtered_step)

   return list_of_lists_of_single_steps

def map_genome_to_modules(ncbi_taxonomy_id):

   ncbi_taxonomy_id = str(ncbi_taxonomy_id)

   module_definitions = open("../ref-dbs/module_definitions.tsv")
   kegg_genome_file   = glob.glob("../ref-dbs/kegg_genomes/" + ncbi_taxonomy_id + "/*.json") 
   kegg_genome_file   = open(kegg_genome_file[0])
   kegg_genome        = json.load(kegg_genome_file)


def break_down_complex_step(step, defn, md):

   if step[0] != "(" or step[-1] != ")": 
      step = "(" + step + ")"

   openings = step.split("(")

   if len(openings) == 2 and openings[0] == "": 

      alternatives = openings[1][:-1]
      return alternatives


   else:


      alternatives = []

      nodes   = []
      ko_term = ''

      open_parenth_counter  = 0
      closed_parent_counter = 0
      priority_starts       = False # referring to the inner ()
      priority_stops        = False # referring to the inner ()
      alternative_option    = False
      semi_step             = False

      tailed  = step[1:-1]
      copy    = tailed + "*"

      indices_for_unique_alternatives = [0]

      for index, character in enumerate(copy): 

         """
         Loop to get the indices of the unique alternative roads
         of a step.
         """         

         if character == "*":
            continue

         if character == "(": 

            priority_starts = True
            open_parenth_counter += 1
            ko_term = ''
         
         elif character == ";": 

            semi_step = True
            ko_term   = ''

         elif character == ",": 
            ko_term = ''
            alternative_option = True

         elif character == ")": 

            closed_parent_counter += 1
            ko_term = ''
            priority_stops = True
            alternative_option = False

            
         ko_term = character
         node      = ko_term ; nodes.append(node)
         level     = open_parenth_counter - closed_parent_counter
         link_from = copy[copy.index(node)-1]
         link_to   = copy[index + 1]

                  

         if (link_to == "," or link_from == ",") and level <= 0: 

            if len(indices_for_unique_alternatives) > 1: 
               
               if index - indices_for_unique_alternatives[-1] > 3:
                  indices_for_unique_alternatives.append(index + 1)
            
            else:

               indices_for_unique_alternatives.append(index + 1)


      unique_alternatives = split_stirng_based_on_indeces(copy, indices_for_unique_alternatives)
      
      counto = 0
      for alternative in unique_alternatives:
         counto += 1
         # print(counto, alternative)

         if alternative.count("K") == 1:
            alternative = alternative.replace(",", "")
            alternative = alternative.replace("*", "")
            alternative = alternative.replace(";", "")
            
            alternatives.append(alternative)


         else:

            print(counto, alternative)

         # if copy[index + 1] in (";", "(", ")", ",", "*"):


         #    print("link_from: ", link_from)
         #    print("link_to:", link_to)

         #    print("priority_starts: ", priority_starts)
         #    print("priority_stops: ", priority_stops)
         #    print("alternative_option: ", alternative_option)
         #    print("semi_step: ", semi_step)
         #    print("level: ", level)

         #    print("\n---\n")



            # if level == 0 and alternative_option == True: 

            #    alternatives.append([node])


            # if len(alternatives) == 0 :
            #    alternatives.append([node])

            # if semi_step: 

            #    for alternative in alternatives:

            #       alternative.append(node)
                  
            #       semi_step = False
            #       alternative_option = False


            # if alternative_option and level == 0: 
            #    alternatives.append([node])


            # if alternative_option and priority_starts and priority_stops == False: 

            #    alternatives.append([node])
            #    alternative_option = False 



         
      print(copy, md)
      print("\n~~~~\n")
      # sys.exit(0)







def split_stirng_based_on_indeces(s, indices):
   """
   REMEMBER! You need to give '0' as your first index, e.g.:
   indices = [0,5,12,17]
   """
   parts = [s[i:j] for i,j in zip(indices, indices[1:]+[None])]
   return parts


















parse_definitions_file()
# map_genome_to_modules(36033)
