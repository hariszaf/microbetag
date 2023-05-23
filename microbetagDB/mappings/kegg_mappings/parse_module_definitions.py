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
import re, sys
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

def create_final_regular_dictionary(module_steps_parsed, module_components_raw):

    final_regular_dict = {}
    counter = 0

    # Parse module steps and export them into a text file
    for module, steps in module_steps_parsed.items():

        final_regular_dict[module] = {}
        step_number = 0

        # Deal with one step at a time
        for step in steps:

            check = False

            # Build temp_string ----    First part of the function
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
                    if count >= 2 :
                        temp_string += char

                elif char == ")":
                    count -= 1
                    if count >= 1 or (options == 1 and count >=1):
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

            # Steps, i.e. temp_string, with more than one alternatives ----    Second part of the function
            if options >= 2:

                counter += 1

                # Check for "-" terms
                temp_string = re.sub("-%(.*)", "", temp_string, count=0, flags=0)
                if "-" in temp_string:
                    temp_string = re.sub(r'-K[0-9]{5}', '', temp_string)

                # 
                alts_for_a_step = [] ; check = True

                # I am splitting using the space, as it denotes an INDEPENDENT alternative for the current step
                alternatives = temp_string.split()

                # and exam each alternatve separately
                for alt in alternatives:

                    # Use regular expressions to extract the groups within parentheses
                    """
                    it returns the inner parentheses, for example:
                    alt: ((K03831,K03638)_K03750)
                    will return:
                    K03831,K03638
                    """
                    pattern = r"\(([^()]+)\)"
                    matches = re.findall(pattern, alt)

                    alt_ops = []

                    if alt[:2] == "((":
                        alt = alt[1:-1]

                    pts = re.split(r"\)\+|\)_\(|_\(|\)_", alt)

                    for index, pt in enumerate(pts):
                        if ")" in pt or "(" in pt:
                            i = re.sub("\(|\)", "", pt)
                            pts[index] = i

                    for index, pt in enumerate(pts):
                        if "," in pt:
                            d = []
                            for i in pt.split(","):
                                d.append([i])
                            pts[index] = d
                        elif "_" in pt:
                            h = []
                            for j in pt.split("_"):
                                h.append(j)
                            pts[index] = h
                        elif "+" in pt:
                            d = []
                            for i in pt.split("+"):
                                d.append([i])
                            pts[index] = d
                        else:
                            pts[index] = [pt]

                    for index_1, i in enumerate(pts):
                        if isinstance(i,list):
                            for index_2, j in enumerate(i):
                                if isinstance(j,list):
                                    for index_3, z in enumerate(j):
                                        if "+" in z or "_" in z:
                                            pts[index_1][index_2][index_3] = re.split(r"\+|_", z)
                                            pts[index_1][index_2] = pts[index_1][index_2][index_3]

                    alts_for_a_step.append(pts)


            # Single case steps. Most cases are single KOs or "choose one from" (e.g. (K00969,K06210)) or just a combination of KOs (e.g., K06215+K08681)
            else:
                # Check for "-" terms
                if "-" in step:
                    temp_string = [re.sub(r'-K[0-9]{5}', '', step)]

                # Split if ","; if not one-long list is made
                else:
                    temp_string = step.split(",")
                temp_string = [re.sub(r'\(|\)', '', i) for i in temp_string]


                for i in range(len(temp_string)):
                    temp_string[i] = re.split(r"\+|_", temp_string[i])


            # If temp_string is a string, split in spaces
            if isinstance(temp_string, str):
                temp_string = temp_string.split()


            # In all cases, temp_string is a list now. 
            temp_string = sorted(temp_string, key=len)


            # Check if a list is there
            if check: 
                combos = []
                for i in alts_for_a_step:
                    if len(i) > 1:
                        j = list(itertools.product(*i))
                        for case in j:
                            combos.append(list(flatten(case)))
                    else:
                        combos.append(i)
                final_regular_dict[module][step_number] = combos
            else:
                final_regular_dict[module][step_number] = temp_string

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

P = create_final_regular_dictionary(module_steps_parsed, module_components_raw)

q = {}
for md, steps in P.items():

    q[md] = {}
    q[md]["definition"] = module_components_raw[md]
    q[md]["#-of-steps"] = len(steps)    
    q[md]["steps"] = steps


    # for step, cases in steps.items():
    #     if len(q[md]["steps"]) > 0:
    #         q[md]["steps"].append([])
    #     q[md]["steps"] = q[md]["steps"] + []
    #     for case in cases:
    #         alts = re.split(r"\+|_", case)
    #         if len(q[md]["steps"]) == 0:
    #             q[md]["steps"].append([alts])
    #         else:
    #             q[md]["steps"][-1] = q[md]["steps"][-1] + [alts]
    # # We make the set again a list, in order to be able to dump the dictionary to a json file
    # q[md]["unique-KOs"] = [list(set(list(flatten(q[md]["steps"]))))]
    # q[md]["all-alts"] = list(itertools.product(*q[md]["steps"]))






import json
with open('result.json', 'w') as fp:
    json.dump(P, fp, indent=4)


print(P["M00855"]) 
print(q["M00022"])
print(list(flatten(q["M00022"]["steps"]["1"])))


