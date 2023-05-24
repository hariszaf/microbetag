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

def parenthetic_contents(string):
    """Generate parenthesized contents in string as pairs (level, contents)."""
    stack = []
    for i, c in enumerate(string):
        if c == '(':
            stack.append(i)
        elif c == ')' and stack:
            start = stack.pop()
            yield (len(stack), string[start + 1: i])

def get_strong_or(string):

    # Build temp_string ----    First part of the function
    count = 0
    options = 0
    parsing_string = ""
    for char in string:

        if char == "(":
            count += 1
            options += 1
            if len(parsing_string) > 1 and parsing_string[-1] == "-":
                parsing_string += "%"
            if count >= 2 :
                parsing_string += char

        elif char == ")":
            count -= 1
            if count >= 1 or (options == 1 and count >=1):
                parsing_string += char
            else:
                continue

        elif char == ",":
            if count >= 2:
                parsing_string += char
            else:
                parsing_string += " "
        else:
            parsing_string += char

    return parsing_string, options

def split_com_and_uper(alist):

    try:
        for i in range(len(alist)):

            if "," in alist[i]:

                alist[i] = alist[i].split(",")

        for i in range(len(alist)):

            for j in range(len(alist[i])):

                if "+" in alist[i][j]:
                    alist[i][j] = alist[i][j].split("+")
    except:
        pass

    return alist

def disentangle_splited_alternative(list_of_pars, level, true_alt, md):

    check = False
    while level > 0:
        for par1 in list_of_pars: 
            if par1[0] == level: 
                for par2 in list_of_pars: 
                    if par2[0] == level - 1 :
                        if par1[1] in par2[1]:
                            check = True

                            index = par2[1].index(par1[1])

                            # [IMPORTANT] NOW I KNOW THAT BOTH par1 and par2 ARE INNER TO THE true_alt PARENTHESIS
                            if par2[0] > 0:
                                print(par2, "`````", par1, "~~~~~~~", true_alt)


                                continue
                            # THE REST, YOU TREAT LIK ANY OTHER PART
                            else:
                                print(">>", true_alt)

                            # print(true_alt, "~~", get_strong_or(true_alt[0]), "~",  md)


        level -= 1

    if not check:
        if "_" not in true_alt:

            c = re.split(r"\)\+|\+\(", true_alt)
            d = [s.strip('(|)') for s in c] 
            e = split_com_and_uper(d)

            if any(isinstance(i, list) for i in e):
                e = [[i] if isinstance(i,str) else i for i in e  ]
                e = list(itertools.product(*e))
                return (list(itertools.product(*e)))
            else: 
                e = [e[0].split("+")]
                return e
        else:
            c = (true_alt.split("_")) ; d = []
            for i in range(len(c)):
                if c[i][0]=="(" and "," in c[i] and ")" not in c[i]:
                    d.append(c[i] + "_"+ c[i+1])
                else:
                    d.append(c[i])
            d = [s.strip('(|)') for s in d]
            e = split_com_and_uper(d)
            if any(isinstance(i, list) for i in e):
                e = [[i] if isinstance(i,str) else i for i in e  ]
                e = list(itertools.product(*e))
                # print(e)
                return e
            else: 
                return e





def parse_regular_module_dictionary(module_components_raw, structural_list):
    """
    bifurcating_list: consider adding it or not, if yes give it as an argument
    """
    # Parse raw module information
    module_steps_parsed = {}
    for key, values in module_components_raw.items():
        # if key != "M00022":
        #     continue
        values = values.replace(" --", "")
        values = values.replace("-- ", "")
        
        # or key in bifurcating_list 
        if key in structural_list:
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

def create_final_regular_dictionary(module_steps_parsed):
    """
    This function returns all the possible combinations of KOs to have a complete KEGG module
    """

    final_regular_dict = {}

    # Parse module steps and export them into a text file
    for module, steps in module_steps_parsed.items():

        final_regular_dict[module] = {}
        step_number = 0 

        # Deal with one step at a time
        for step in steps:

            step_number += 1
            check = False

            temp_string, options = get_strong_or(step)

            
            # Steps, i.e. temp_string, with more than one alternatives ----    Second part of the function
            if options > 1:

                # Check for "-" terms
                if "-%" in temp_string:
                    temp_string = re.sub("-%\(.*\)", "", temp_string, count=0, flags=0)
                if "-" in temp_string:
                    temp_string = re.sub(r'-K[0-9]{5}', '', temp_string)

                # I am splitting using the space, as it denotes an INDEPENDENT alternative for the current step
                alternatives = temp_string.split()

                # Keep all alternatives of a step in a list
                alts_for_a_step = [] ; check = True

                # Parse each alternative to get the possible combinations of it
                for alternative in alternatives:

                    splited_alternative = list(parenthetic_contents(alternative))

                    # add the alternative as a list of KO terms in the alternatives list
                    if len(splited_alternative) == 0:
                        """
                        the alternative on the right becomes like in the left
                        [['K00665']] ~~~ K00665
                        [['K00665'], ['K00667', 'K00668']] ~~~ K00667+K00668
                        """
                        sp = re.split(r"\+|_", alternative)
                        alts_for_a_step.append(sp)

                    else:
                        max_levels = max([o[0] for o in splited_alternative])

                        q = disentangle_splited_alternative(splited_alternative, max_levels, alternative, module)








                # # and exam each alternatve separately
                # for alt in alternatives:

                #     if alt[:2] == "((":
                #         alt = alt[1:-1]

                #     # pts = re.split(r"\)_\(", alt)

                #     pts = re.split(r"\)\+|\)_\(|_\(|\)_|_\(\(", alt)


                #     for index, pt in enumerate(pts):
                #         if ")" in pt or "(" in pt:
                #             i = re.sub("\(|\)", "", pt)
                #             pts[index] = i                        
                            

                #         if "," in pt:
                #             d = []
                #             for i in pt.split(","):
                #                 d.append([i])
                #             pts[index] = d

                #         elif "_" in pt:
                #             h = []
                #             for j in pt.split("_"):
                #                 h.append(j)
                #             pts[index] = h


                #         elif "+" in pt:
                #             d = []
                #             for i in pt.split("+"):
                #                 d.append([i])
                #             pts[index] = d
                #         else:
                #             pts[index] = [pt]


                #     for index_1, i in enumerate(pts):
                #         if isinstance(i,list):
                #             for index_2, j in enumerate(i):
                #                 if isinstance(j,list):
                #                     for index_3, z in enumerate(j):
                #                         if "+" in z or "_" in z:
                #                             pts[index_1][index_2][index_3] = re.split(r"\+|_", z)
                #                             pts[index_1][index_2] = pts[index_1][index_2][index_3]

                #     alts_for_a_step.append(pts)

            # Single case steps. Most cases are single KOs or "choose one from" (e.g. (K00969,K06210)) or just a combination of KOs (e.g., K06215+K08681)
            else:

                # Check for "-" terms
                if "-" in step:
                    # print(temp_string, "~~~", step)
                    temp_string = [re.sub(r'-K[0-9]{5}', '', step)]
                    # print(temp_string, "~~~", step)

                # Split if ","; if not one-long list is made
                else:
                    temp_string = step.split(",")
                temp_string = [re.sub(r'\(|\)', '', i) for i in temp_string]


                # [ALERT] This must lead to errors... 
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


module_steps_parsed = parse_regular_module_dictionary(module_components_raw, structurals)

P = create_final_regular_dictionary(module_steps_parsed)

q = {}
for md, steps in P.items():

    q[md] = {}
    q[md]["id"] = md
    q[md]["definition"] = module_components_raw[md]
    q[md]["#-of-steps"] = len(steps)    
    q[md]["steps"] = steps
    q[md]["unique-KOs"] = list(set(list(flatten(q[md]["steps"].values()))))

    # Make sure all alternatives of a step have the same format
    for step, opt in q[md]["steps"].items():
        for index, c in enumerate(opt):
            if len(c) == 1:

                if len(c[0]) == 1:
                    
                    flat_list = [item for sublist in c for item in sublist]
                    q[md]["steps"][step][index] = flat_list 

                else:
                    x = all(len(i) == 1 for i in c[0])
                    if x:

                        flat_list = [item for sublist in c for item in sublist]
                        if flat_list:
                            if flat_list[0] != "K":
                                q[md]["steps"][step][index] = flat_list

for md, steps in q.items():
    for step, opts in q[md]["steps"].items():
        for index, opt in enumerate(opts):
            if isinstance(opt[0], list):
                q[md]["steps"][step][index] = list(flatten(opt))


import json
# with open('result.json', 'w') as fp:
#     json.dump(q, fp, indent=4)


