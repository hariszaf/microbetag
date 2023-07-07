import sys, os
import json, re
from pathlib import Path
import itertools

def write_sequential_ko_pairs(sequence, outfile):

    for ko1, ko2 in zip(sequence, sequence[1:]):
        outfile.write(ko1 + "\t" + ko2 + "\n")

    return 



def overall_module_kos_graph(mos_kos_map):
    """
    This function builds/writes a 2-column file based on the module_definition_map.json 
    under mappings/kegg_mappings, whith all the paired-end steps described on that, 
    meaning for example,
    M00093 with the definition : K00981 (K00998,K17103) K01613 will return:
    
    K00981  K00998
    K00981  K17103
    K00998  K01613
    K17103  K01613
    """

    kos_numb_of_mod = {}

    mos_kos_graph = open("kos_on_modules_graph.tsv", "w")

    for module, data in mos_kos_map.items():

        if module == "md:M00009":

            for step_a, step_b in zip(data["steps"], data["steps"][1:]):


                print("step_a: ", step_a)

                print("step_b: ", step_b)

                if isinstance(step_a, str): 
                    step_a = [step_a]
                if isinstance(step_b, str):
                    step_b = [step_b]

                types_step_a = [type(item) for item in step_a]
                types_step_b = [type(item) for item in step_b]

                if len(types_step_a) == 1 and types_step_a[0] == str and len(types_step_b) == 1 and types_step_b[0] == str: 
                    print("hello friend")
                    mos_kos_graph.write(step_a[0] + "\t" + step_b[0] + "\n")

                else:
                    fromKO = set()
                    for case in range(len(types_step_a)):
                        if types_step_a[case] == list:
                            write_sequential_ko_pairs(step_a[case], mos_kos_graph)
                            fromKO.add(step_a[case][-1])
                        else:
                            fromKO.add(step_a[case])

                    toKO = set()
                    for case in range(len(types_step_b)):
                        if types_step_b[case] == list:
                            print("--> ", step_b[case])
                            write_sequential_ko_pairs(step_b[case], mos_kos_graph)
                            toKO.add(step_b[case][0])
                        else:
                            toKO.add(step_b[case])

                    for pair in list(itertools.product(fromKO, toKO)):
                        mos_kos_graph.write(pair[0] + "\t" + pair[1] + "\n")


    return 1



def build_kos_graph(genome_kos):

    """
    Aim of this function is to build the graph of KOs of a genome related to modules,
    in the format NetCooperate can use it as input, i.e.:
    K00001  K00002
    K00002  KO0003
    
    This function will be performed to all genomes under ref-dbs/all_genomes
    """


    return genome_graph





def parse_module_definition(definitions):

    module_steps_parsed = {}
    for module_id, definition in definitions.items():
        definition = definition.replace(" --", "")
        definition = definition.replace("-- ", "")

        module = []
        parenthesis_count = 0
        for character in definition:
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
        module_steps_parsed[module_id] = steps

        # Remove modules depending on other modules
        temporal_dictionary = module_steps_parsed.copy()
        for key, values in temporal_dictionary.items():
            for value in values:
                if re.search(r'M[0-9]{5}', value) is not None:
                    del module_steps_parsed[key]
                    break
        print(module)
    return module_steps_parsed









if __name__ == "__main__":
    upath = Path(os.getcwd())
    mappings_file = upath.parent.absolute() / "mappings/kegg_mappings/module_definitions.tsv"  #module_definition_map.json"
    f = open(mappings_file, "r")
    
    definitions_dic = {}
    for line in f:
        module, definition = line.split("\t")[0], line.split("\t")[1]
        definitions_dic[module] = definition

    parsed_modules = parse_module_definition(definitions_dic)

    for i, j in parsed_modules.items():
        print(i)
        print(definitions_dic[i])
        print(j)
        print("~~~")


    sys.exit(0)
    
    gmap = json.load(f)
    overall_module_kos_graph(gmap)





