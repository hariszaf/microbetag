#!/usr/bin/env python3

"""
Aim of this script is to build all the unique sets 
of KO terms that can be used to build up each KEGG 
module (https://www.genome.jp/brite/ko00002)

For example: 
definition of md:M00545: ((K05708+K05709+K05710+K00529);K05711),K05712
would return 2 sets:
1. K05708, K05709, K05710, K00529, K05711
2. K05712

while 
definition of md:M00546:  K13484,K07127;(K13485,K16838,K16840)
would return: 
1. K13484
2. K07127, K13485
3. K07127, K16838
4. K07127, K16838

The json file will have 2 levels, in the first one, the various steps will
be denoted and in the second the multiple alternative combinations of terms 
will be shown. 
All the terms of a combination will be needed for the module to be complete
"""

import glob
import json
import sys
import itertools
import re
from collections.abc import Iterable


def split_definition_to_steps(definition, md):

   """
   Split where there is a ';' being outside of a parenthesis
   """

   quasi_steps           = definition.split(";")
   count_front_parenth   = 0
   count_reverse_parenth = 0
   pseudo_step           = []
   in_parenthesis        = False
   complete_steps        = []

   for quasi_step in quasi_steps:

      if "(" not in quasi_step and ")" not in quasi_step and in_parenthesis == False:

            modules = []

            if "," in quasi_step: 

               quasi_step = quasi_step.split(",")
               for x in quasi_step:
                  modules.append(x)

            else:
               modules = [quasi_step]

            if modules:
               complete_steps.append(modules)

      else:

         if "(" in quasi_step or ")" in quasi_step:

            num_of_open            = quasi_step.count('(')
            count_front_parenth   += num_of_open
            num_of_closed          = quasi_step.count(')')
            count_reverse_parenth += num_of_closed

         if count_front_parenth > count_reverse_parenth:
            in_parenthesis = True
            pseudo_step.append(quasi_step)
                     
         elif count_front_parenth == count_reverse_parenth and in_parenthesis == False:

            complete_steps.append(quasi_step)

         elif count_front_parenth == count_reverse_parenth and in_parenthesis == True:
            in_parenthesis = False
            pseudo_step.append(quasi_step)
            complete_step = ' '.join(pseudo_step)
            complete_steps.append(complete_step)
            pseudo_step = []

   if complete_steps:
      return complete_steps

def parse_definitions_file():

   definitions_file = open("../ref-dbs/module_definitions.tsv", "r")
   module_definitions_steps = {}
   

   for line in definitions_file:

      md           = line.split("\t")[0]
      definition   = line.split("\t")[1][:-1]
      parsed_steps = []

      # if md != "md:M00785_2":
      #    continue

      if "(" not in definition and "," not in definition:

         """
         no alternatives per step 
         REMEMBER! A step can be a complex, so more than 1 KOs may be required for a single sptep, 
         meaning that a single step (a list) may include multiple lists 
         """

         if ";" in definition:
            steps = definition.split(";")
            parsed_steps = parse_valid_steps_of_a_module(steps)

         if ";" not in definition:
            if "+" in definition:
               parsed_steps = [[definition.split("+")]]

         if ";" not in definition: 
            parsed_steps = [[definition]]

      else:

         complete_steps = split_definition_to_steps(definition, md)

         for step in complete_steps:

            if isinstance(step, list):
               parsed_steps.append(step)
            
            else:
               alternatives_for_a_step = break_down_complex_step(step, definition, md)
               parsed_steps.append(alternatives_for_a_step)

         pseudo_step           = []
         count_front_parenth   = 0
         count_reverse_parenth = 0
         in_parenthesis        = False


      module_definitions_steps[md]               = {}
      y                                          = 0
      num_of_steps                               = len([y+1 for x in parsed_steps if "--" not in x ])
      module_definitions_steps[md]['#-of-steps'] = num_of_steps

      parsed_steps                               = tuple(parsed_steps)
      module_definitions_steps[md]['steps']      = parsed_steps

      unique_ko_terms = set(list(flatten(parsed_steps)))
      module_definitions_steps[md]['unique-KOs'] = list(unique_ko_terms)


   """
   Here we save our actual output as a .json file
   """

   with open("../ref-dbs/module_definition_map.json", "w") as f:
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
   
         semis = []
         step  = step.split("-")
         semi_counter = 0
         
         if len(step) == 2 and step[0] != "":
            if "+" in step[0]:
               components = step[0].split("+")
               list_of_lists_of_single_steps.append(components)

         else: 
            # Step has more than 2 parts
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

            list_of_lists_of_single_steps.append([filtered_step])  # check if needed the list!! 

   return list_of_lists_of_single_steps

def break_down_complex_step(step, defn, md):

   if step[0] != "(" or step[-1] != ")": 
      step = "(" + step + ")"

   if "-" in step:
      minus_indices = [0,]
      [minus_indices.append(i) for i, c in enumerate(step) if c == "-"]


      step_with_no_minus = split_stirng_based_on_indeces(step, minus_indices)
      parts_of_no_minus_to_keep = [step_with_no_minus[0]]


      for entry in step_with_no_minus[1:]: 
         # ORIGINAL ! 
         # check = [e for e in [")", "(", "+", ";"] if e in entry]

         check = [entry.index(e) for e in [")", "(", "+", ";"] if e in entry]

         if len(check) == 0:
            continue
         else:
            # ORIGINAL - QUITE WORKING!!! 
            # continue_from = entry.index(check[0])
            continue_from = sorted(check)[0]
            print("the entry: ", entry)
            print("check: ", check)
            print("continue: ", continue_from)
            parts_of_no_minus_to_keep.append(entry[continue_from:])

      step = ''.join(parts_of_no_minus_to_keep)
      print("STEP::: ", step)

   openings = step.split("(")

   if len(openings) == 2 and openings[0] == "": 

      alternatives = []
      options = openings[1][:-1]

      if "," in options:
         alternatives = options.split(",")

         for case in range(len(alternatives)):
            if " " in alternatives[case]:
               alternatives[case] = alternatives[case].split(" ")

      else:
         alternatives.append(options)

      alternatives = [x for x in alternatives if x]
      print([x for x in alternatives if x])
      return alternatives


   else:

      alternatives          = []
      node                  = ''
      ko_term               = ''
      open_parenth_counter  = 0
      closed_parent_counter = 0
      node_index            = -1
      alternative_option    = False
      semi_step             = False

      for character in step[:-1]:
         if character == "(":
            open_parenth_counter += 1
         if character == ")":
            closed_parent_counter += 1

      if open_parenth_counter > closed_parent_counter:
         """
         then the last character of copy needs to be ')'
         """
         step  = step[1:-1]

      copy    = step + "*"

      # Here is the list where we append indices of the string 
      # under study to split it the way it fits us
      indices_for_unique_alternatives = [0]
      open_parenth_counter  = 0
      closed_parent_counter = 0

      print("COPY: ", copy)

      for index, character in enumerate(copy): 

         """
         Loop to get the indices of the unique alternative roads
         of a step.
         """

         if character == "*":
            continue

         (
            open_parenth_counter,
            closed_parent_counter,
            semi_step,
            alternative_option, 
            ko_term,
            link_from,
            link_to,
            level,
            node_index, node
         )                    = handle_character(
                                                   copy,
                                                   index,
                                                   character, 
                                                   open_parenth_counter, 
                                                   closed_parent_counter, 
                                                   semi_step, 
                                                   alternative_option, 
                                                   ko_term,
                                                   node_index,
                                                   node
                                                )

         if (link_to == "," or link_from == ",") and level <= 0: 

            """
            Level denotes how many parentheses have been opened and closed.
            If their diffrence equals to zero, then it means we have distinct alternatives
            """

            if len(indices_for_unique_alternatives) > 1: 
               
               if index - indices_for_unique_alternatives[-1] > 3:
                  indices_for_unique_alternatives.append(index + 1)

            else:

               indices_for_unique_alternatives.append(index + 1)
      
      node_index          = -1
      unique_alternatives = split_stirng_based_on_indeces(copy, indices_for_unique_alternatives)

      open_parenth_counter  = 0 
      closed_parent_counter = 0

      print("UNIQUE: ", unique_alternatives)

      for index, alternative in enumerate(unique_alternatives):

         alternative = alternative.replace("*", "")

         if alternative.count("K") == 1:

            alternative = alternative.replace(",", "")
            alternative = alternative.replace(";", "")
            alternatives.append(alternative)

         else:

            if alternative and alternative[0] == ",": 
               alternative = alternative[1:]

            if "," not in alternative:

               alternative = alternative.replace("(", "_")
               alternative = alternative.replace(")", "_")
               alternative = alternative.replace("+", "_")
               alternative = alternative.replace(";", "_")
               alternative = alternative.split("_")               
               alternative = [i for i in alternative if i]

               alternatives.append(alternative)

            else:
               
               opening_parenth = 0
               closing_parenth = 0
               complet_parenth = 0

               # remove "()" from start and end if not necessary
               print("SOONER  ", alternative)
               for character in alternative:
                  if character == "(":
                     opening_parenth += 1
                  if character == ")":
                     closing_parenth += 1

                  level = opening_parenth - closing_parenth

                  if level == 0:
                     complet_parenth += 1 

               if complet_parenth == 1:
                  alternative = alternative[1:-1]
                  
               indices  = [i for i, c in enumerate(alternative) if c == ";"]
               openings = [i for i, c in enumerate(alternative) if c == "("]
               closings = [i for i, c in enumerate(alternative) if c == ")"]
               open_close_pairs = [[x,y] for x,y in zip(openings, closings)]
               
               tmp          = []
               break_points = [0,]
               counter = 0 
               for index in indices:
                  break_point = True
                  for pair in open_close_pairs: 
                     if index > pair[0] and index < pair[1]:
                        break_point = False

                  if break_point: 
                     break_points.append(index)

               print(">>>> alternative: ", alternative, md)
               clean = parse_to_a_list_of_single_items(alternative)
               print(">>>> Clean: ", clean)
               pools = get_possible_alternatives_from_a_step(clean)
               print(">>>> POols: ", pools)

               alternatives = combine_alternatives(alternatives, pools)


               print("!!!!! Alternatives : ", alternatives)
               alternatives = [list(flatten(i)) for i in alternatives if isinstance(i, list)]
               print("!!!!! Flattened lternatives : ", alternatives)



               print(alternatives)

      return alternatives

def count_nested_lists(l):
   count = 0 
   for e in l: 
      if isinstance(e, list):
         count = count + 1 + count_nested_lists(e)
   return count

def get_possible_alternatives_from_a_step(clean_step):

   pools = []
   add_extra_paths = False
   complex_present = False

   for index, entry in enumerate(clean_step):

      if entry == " ":
         continue

      if entry == "(":
         if add_extra_paths == False:
            """
            This if statement is crucial to keep track of what is happening inside an 
            already open parenthesis.
            """
            pool = []
            add_extra_paths = True
            continue

      if entry == "+":
         if add_extra_paths:
            complex_present = True
            continue
         else:
            
            continue
      
      if entry == ")":
         if pool:
            pools.append(pool)
            add_extra_paths = False
            pool = []
         continue

      if add_extra_paths == False and complex_present == False and "K" in entry:
         pools.append([entry])

      elif add_extra_paths and complex_present == False and "K" in entry: 
         pool.append(entry)

      elif add_extra_paths == False and complex_present and "K" in entry:
         pools.append([entry])

      elif add_extra_paths and complex_present and "K" in entry:
         if pool: 
            pool[-1] = pool[-1] + "+" + entry

         else:
            pool.append([entry])
         complex_present = False
      
   return pools

def combine_alternatives(alternatives, combos):

   all_combinations = []

   if len(combos) == 1: 
      all_combinations = [tuple(combos[0])]

   elif len(combos) == 2:
      all_combinations = list(itertools.product(combos[0], combos[1]))

   else:
      all_combinations = list(itertools.product(combos[0], combos[1]))
      for lista in combos[2:]:
         all_combinations = list(itertools.product(all_combinations, lista))

   for i in all_combinations:

      single_combo = []
      i = list(i)

      for y in i: 
         
         if isinstance(y, str):
            single_combo.append(y)
         
         else:

            y = list(y)
            for k in y: 

               if isinstance(k, str):
                  single_combo.append(k)

               elif isinstance(k, tuple):
                  single_combo.append(list(flatten(k)))

      if single_combo:
         # remove_pluses = []
         # for j in single_combo:
         #    if "+" in j:
         #       remove_pluses.append(j.split("+"))
         #    elif " " in j:
         #       remove_pluses.append(j.split(" "))
         #    else:
         #       remove_pluses.append(j)
         # remove_pluses = list(flatten(remove_pluses))
         alternatives.append(single_combo)
         # alternatives.append(remove_pluses)
   
   return alternatives

def flatten(lis):
     for item in lis:
         if isinstance(item, Iterable) and not isinstance(item, str):
             for x in flatten(item):
                 yield x
         else:        
             yield item

def parse_to_a_list_of_single_items(rule):

   parts_of_alts = rule.split("K")
   parts_of_rule = []
   
   for i in parts_of_alts: 

      if any(map(str.isdigit, i)):
         new_i = "K"+i
         asdas = re.split('([^a-zA-Z0-9])', new_i)
         parts_of_rule.append([x for x in asdas if x != ""])
         
      elif i == "":
         continue

      else:
         parts_of_rule.append(i)

   clean = []

   for x in parts_of_rule:
      if isinstance(x,str):
         clean.append(x)
      else:
         for y in x:
            clean.append(y)
   
   # with respect to the M00854 case
   if clean[0] == "(" and (clean.count("(") != clean.count(")")):
      clean = clean[1:]

   return clean

def handle_character(part_of_def, index, character, open_parenth_counter, closed_parent_counter, 
                     semi_step, alternative_option, ko_term, node_index, node):

   if character == "(": 
      open_parenth_counter += 1
      ko_term = ''

   elif character == ";": 
      semi_step = True
      ko_term   = ''

   elif character == ",": 
      alternative_option = True
      ko_term = ''

   elif character == ")": 
      closed_parent_counter += 1
      ko_term = ''
      alternative_option = False

   elif character == "K":
      node_index = index

   level     = open_parenth_counter - closed_parent_counter
   link_to   = part_of_def[index + 1]


   if index == 0 and node_index < 0:
      link_from = '' ; 

   elif index > 0 and node_index < 0: 
      link_from = part_of_def[index - 1]

   elif node_index == 0:
      """
      part_of_def starting from K0, e.g.: K01196,((K00705,K22451);(K02438,K01200))
      """
      link_from = "^"

   else:
      link_from = part_of_def[index - 1]

   if link_from and link_to in (";", "(", ")", ",", "^"):
      node = part_of_def[node_index:index + 1]

   return open_parenth_counter, closed_parent_counter, semi_step, alternative_option,\
            ko_term, link_from, link_to, level, node_index, node

def split_stirng_based_on_indeces(s, indices):
   """
   REMEMBER! You need to give '0' as your first index, e.g.:
   indices = [0,5,12,17]
   """
   parts = [s[i:j] for i,j in zip(indices, indices[1:]+[None])]
   return parts

parse_definitions_file()

