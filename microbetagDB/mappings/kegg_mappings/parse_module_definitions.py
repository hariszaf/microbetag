"""
The pathway module is defined by the logical expression of K numbers, and the signature module 
is defined by the logical expression of K numbers and M numbers.
A SPACE ( ) or a PLUS (+) sign, representing a connection in the pathway or the molecular complex, is treated as an AND operator 
and a COMMA (,), used for alternatives, is treated as an OR operator. 
A MINUS (-) sign designates an optional item in the complex. 

Aim of this script is to build all the unique sets 
of KO terms that can be used to build up each KEGG 
module (https://www.genome.jp/brite/ko00002)


This script is based on 2 function from the following script of microbeAnnotator:
https://github.com/cruizperez/MicrobeAnnotator/tree/master/microbeannotator/data/01.KEGG_DB/00.KEGG_Data_Scrapper.py



The json file will have 2 levels, in the first one, the various steps will
be denoted and in the second the multiple alternative combinations of terms 
will be shown. 
All the terms of a combination will be needed for the module to be complete
"""

import itertools
import re
import collections.abc
collections.Iterable = collections.abc.Iterable



bifs = [
         "M00373",
         "M00532",
         "M00376",
         "M00378",
         "M00088",
         "M00031",
         "M00763",
         "M00133",
         "M00075",
         "M00872",
         "M00125",
         "M00119",
         "M00122",
         "M00827",
         "M00828",
         "M00832",
         "M00833",
         "M00837",
         "M00838",
         "M00785",
         "M00307",
         "M00048",
         "M00127",
         "M00893",
         "M00895",
         "M00896",
         "M00897",
         "M00898",
         "M00899",
         "M00911",
         "M00913",
         "M00917",
         "M00935"
      ]

structurals = [
   "M00144",
   "M00149",
   "M00151",
   "M00152",
   "M00154",
   "M00155",
   "M00153",
   "M00156",
   "M00158",
   "M00160"
]

def flatten(lis):
     for item in lis:
         if isinstance(item, collections.Iterable) and not isinstance(item, str):
             for x in flatten(item):
                 yield x
         else:        
             yield item

def parse_regular_module_dictionary(bifurcating_list, structural_list, module_components_raw):

    # Parse raw module information
    module_steps_parsed = {}
    for key, values in module_components_raw.items():
        # if key != "M00022":
        #     continue
        values = values.replace(" --", "")
        values = values.replace("-- ", "")
        if key in bifurcating_list or key in structural_list:
            continue
        else:
            module = []
            parenthesis_count = 0
            for character in values:
                if character == "(":
                    parenthesis_count += 1
                    module.append(character)
                elif character == " ":
                    if parenthesis_count == 0:
                        module.append(character)
                    else:
                        module.append("_")
                elif character == ")":
                    parenthesis_count -= 1
                    module.append(character)
                else:
                    module.append(character)
            steps = ''.join(module).split()
            module_steps_parsed[key] = steps

    # Remove modules depending on other modules
    temporal_dictionary = module_steps_parsed.copy()
    for key, values in temporal_dictionary.items():
        for value in values:
            if re.search(r'M[0-9]{5}', value) is not None:

                del module_steps_parsed[key]
                break

    return module_steps_parsed

def create_final_regular_dictionary(module_steps_parsed, module_components_raw, outfile):
    final_regular_dict = {}
    # Parse module steps and export them into a text file
    with open(outfile, 'w') as output:
        for key, value in module_steps_parsed.items():
            output.write("{}\n".format(key))
            output.write("{}\n".format(module_components_raw[key]))
            output.write("{}\n".format(value))
            output.write("{}\n".format("=="))
            final_regular_dict[key] = {}
            step_number = 0

            for step in value:
               step_number += 1
               count = 0
               options = 0
               temp_string = ""
               for char in step:
                  if char == "(":
                        count += 1
                        options += 1
                        if len(temp_string) > 1 and temp_string[-1] == "-":
                           temp_string += "%"
                  elif char == ")":
                        count -= 1
                        if count >= 1:
                           temp_string += char
                        else:
                           continue
                  elif char == ",":
                        if count >= 2:
                           temp_string += char
                        else:
                           temp_string += " "
                  else:
                        temp_string += char
               # print(temp_string, "<<")

               if options >= 2:

                  temp_string = temp_string.replace(")_", "_")

                  if re.search('\)\+', temp_string) is not None:

                     new_tmp = temp_string.split()
                     new_tmp_parts = ""
                     to_remember = []
                     for term in new_tmp: 
                        if re.search('\)\+', term) is not None:
                           parts = term.split(")+")
                        else:
                           to_remember.append(term)
                           continue
                        
                        left = parts[0].split(",")
                        right = parts[1]

                        if isinstance(right, str):
                           right = [right]
                        t = [left, right]
                        combination = [p for p in itertools.product(*t)]
                        for combo in combination:
                           ccombo = ""
                           for element in combo:
                              if len(ccombo)>0:
                                 ccombo =  ccombo + "+" + element
                              else:
                                 ccombo = element
                           new_tmp_parts += ccombo + " "
                     for i in to_remember:
                        new_tmp_parts += i
                     temp_string = new_tmp_parts


                  if re.search('%.*\)', temp_string) is None:
                        temp_string = temp_string.replace(")", "")


                  temp_string = "".join(temp_string.rsplit("__", 1))
                  temp_string = temp_string.split()

               if isinstance(temp_string, str):
                  temp_string = temp_string.split()

               temp_string = sorted(temp_string, key=len)
               final_regular_dict[key][step_number] = temp_string

               output.write("{}\n".format(temp_string))
            output.write("{}\n".format("++++++++++++++++++"))
    return final_regular_dict


modules = open("module_definitions.tsv", "r")
module_components_raw = {}
for line in modules:
   # Remove the "md:" prefix from the id and the new line from the definition
   md, definition = line.split("\t")[0][3:], line.split("\t")[1][:-1]
   # Replace ";" character with a space; this denotes a next layer of the module
   definition = definition.replace(";", " ")
   module_components_raw[md] = definition



module_steps_parsed = parse_regular_module_dictionary(bifs, structurals, module_components_raw)



P = create_final_regular_dictionary(module_steps_parsed, module_components_raw, "MORELES")


# Remove optional KO terms
for md, def_steps in P.items():
      check = False
      for step, combos in def_steps.items():

         for index, case in enumerate(combos):

            if "-%" in case:

               cut_start = case.index("-%")
               indices_object = re.finditer(pattern="\)", string=case)
               indices = [index.start() for index in indices_object]
               new = indices.copy()
               new.append(cut_start)
               new = sorted(new)
               cut_end = new[new.index(cut_start)+1]
               new_case = case[:cut_start] + case[cut_end+1:]
               combos[index] = new_case
               check = True
               continue

            if "-" in case:
               new_case = re.sub(r'-K[0-9]{5}', '', case)
               combos[index] = new_case
               check = True
               continue

         if check: 
            P[md] = def_steps


q = {}

for md, steps in P.items():

   q[md] = {}
   q[md]["#-of-steps"] = len(steps)
   q[md]["steps"] = []

   print("~~~~~~~~~~~~~~")
   print(md)
   print(steps)

   for step, cases in steps.items():

      print(step, cases)

      if len(q[md]["steps"]) > 0:
         q[md]["steps"].append([])


      q[md]["steps"] = q[md]["steps"] + []

      print(">>>>", len(q[md]["steps"]), step)

      for case in cases:

         alts = re.split(r"\+|_", case)


         if len(q[md]["steps"]) == 0:
            q[md]["steps"].append([alts])
         else:
            q[md]["steps"][-1] = q[md]["steps"][-1] + [alts]

   q[md]["unique-KOs"] = [set(list(flatten(q[md]["steps"])))]


print(q["M00011"])


