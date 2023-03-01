#!/usr/bin/env python

import sys, json
import glob
import itertools
# from utils.variables import * 


KMODULES_DEFINITIONS_PARSED = "/home/luna.kuleuven.be/u0156635/github_repos/microbetag/mappings/kegg_mappings/module_definition_map.json"
ALL_GENOMES_MODULES         = "/home/luna.kuleuven.be/u0156635/github_repos/microbetag/ref-dbs/all_genome_modules"

def flatten(list_of_lists):
   """
   This function takes a list of lists and flattens it until it returns a list with 
   all the components of the initial one.
   """
   if len(list_of_lists) == 0:
      return list_of_lists
   if isinstance(list_of_lists[0], list):
      return flatten(list_of_lists[0]) + flatten(list_of_lists[1:])
   return list_of_lists[:1] + flatten(list_of_lists[1:])

def check_for_spaces(l):
   """
   In the modules definition file, there are some cases where complexes are denoted with several 
   KO terms in a single 
   """
   tmp = []
   for term in flatten(l):
      terms = term.split(" ")
      for term in terms: 
         tmp.append(term)
   return tmp

def count_nested_lists(l):
   count = 0 
   for e in l: 
      if isinstance(e, list):
         count = count + 1 + count_nested_lists(e)
   return count

def check_for_minus(definition):
   """ 
   KEGG module definitions may include terms that are not mandatory.
   This function returns a definition list without those.
   """
   for index, list_term in enumerate(definition):
      for case in flatten(list_term):
         if "-" in case:
            definition[index : (index + 2)] = []
   return definition

def split_pluses(path):
   """ 
   In modules with terms in their definitions like 'K02821+K02822+K03475, 
   this functions splits them so it is easier to check whether all the terms are present
   """
   split_path = []
   for entry in path: 
      if "+" in entry:
         entries = entry.split("+")
      
         for i in entries:
            split_path.append(i)
      else:
         split_path.append(entry)

   return split_path

def commonelem_set(z, x):
   one = set(z)
   two = set(x)
   if (one & two):
      print("fuunction") ; print(one & two)
      return True
   else:
      return False # There are no common elements


def extract_alternatives_to_gapfill(ncbi_id_beneficiary, mo_map):
   """
   This function aims at finding KO terms missing for the beneficiary species, 
   that if supplied, could complete a module, based on all the possible alternative paths
   based on the module's definition. 
   """
   list_of_alternatives_to_gap = []
   beneficiarys_genomes = []
   for gfile in glob.glob(ALL_GENOMES_MODULES + "/" + str(ncbi_id_beneficiary) + "/*.json"):
      """
      For a single species/NCBI Taxonomy Id we may have more than one annotated genomes
      [TODO] Consider how we will handle the unity of the alternatives_to_gap
      """
      with open(gfile, 'r') as f:
         genome_mos = json.load(f)

      # Species in-need of by-products
      complete_modules    = []
      alternatives_to_gap = {}

      # Parse genome's KEGG modules and terms as found in the KEGG ORGANISMS db
      counter_complete_modules = 0
      for module, kos_on_its_own in genome_mos.items(): 

         # Get KOs related to the module under stady that are present on the eneficiary's genome
         list_of_kos_present = []
         for k in kos_on_its_own.values():
            if k.startswith("ko:"):
               list_of_kos_present.append(k[3:])
            else:
               list_of_kos_present.append(k)

         # Get the definition of the module under study 
         definition_under_study  = mo_map[module]['steps']
         definition_non_optional = check_for_minus(definition_under_study)

         definition_under_study_proc = []
         for term in definition_non_optional:
            if isinstance(term, str):
               term = [term]
            definition_under_study_proc.append(term)

         """
         The definition_under_study_proc is a list of lists, where any inner list needs to be fully covered
         in order the step to be considered as legit. 
         E.g. [['K00134', 'K00927'], ['K00150', 'K00927'], 'K11389'] 
         means thath we need either both ['K00134', 'K00927'] or both ['K00150', 'K00927'] or just 'K11389'.
         """

         # Get all the possible combinations to have a complete path 
         potential_compl_paths =  [list(tup) for tup in itertools.product(*definition_under_study_proc)]

         # Edit these paths so they are in a flattened format
         flat_potent_compl_paths = []

         for path in potential_compl_paths:
            compl_path = split_pluses(path)
            clear_spaces = check_for_spaces(compl_path)
            flat_potent_compl_paths.append(sorted(clear_spaces))

         # Check whether the kos present on the potential beneficiary species, cover the path
         for path in flat_potent_compl_paths:

            # Check whether all the terms of a path are present in the KOs present on the genome 
            check    =  all(item in list_of_kos_present for item in path)

            # If yes, then the beneficiary does not need anything to perform this module
            if check:
               if module not in complete_modules:
                  counter_complete_modules += 1
                  complete_modules.append(module)

            else:
               # Keep track of the missing steps in the alternative 
               # alteritves_to_gap will look like this:
               #  'md:M00550': {"['K02821', 'K02822', 'K03077', 'K03078', 'K03079', 'K03475', 'K03476']": {'K03079', 'K03476', 'K03078'}}
               # where the values of a module is a dictionary
               # where a path is the key and the missing KOs from the genome is the value

               gaps = set(x for x in set(path) if x not in set(list_of_kos_present))

               # Add the combination of KOs required to be filled to have that module through the path under study
               if module not in alternatives_to_gap:
                  alternatives_to_gap[module] = {}
                  alternatives_to_gap[module][str(path)] = gaps
               else:
                  alternatives_to_gap[module][str(path)] = gaps

      # Remove any complete module
      for key in complete_modules:
         if key in alternatives_to_gap:
            del alternatives_to_gap[key]

      # Keep only the most feasible gap fills 
      for module, val in alternatives_to_gap.items():

         tmp = alternatives_to_gap[module].copy()
         min_val = min([len(val[ele]) for ele in val])
         for path, gaps in alternatives_to_gap[module].items():
            if len(gaps) > min_val + 1:
               del tmp[path]

         alternatives_to_gap[module] = tmp

      list_of_alternatives_to_gap.append(alternatives_to_gap)
      beneficiarys_genomes.append(gfile.split("/")[-1])

   return list_of_alternatives_to_gap, beneficiarys_genomes

def check_for_complements(ncbi_id_donor, list_of_missing_parts, beneficiary_ncbi_id, beneficiayrs_genome_files):
   """ 
   Check whether the potential giver (benefactor), can gap fill one of the alternatives 
   extracted by the extract_alternatives_to_gapfill() function
   """
   counter_donor = 0
   pair_metabolic_interactions = {}
   pair_metabolic_interactions['beneficiary']      = beneficiary_ncbi_id
   pair_metabolic_interactions['donor']            = ncbi_id_donor

   for gfile in glob.glob(ALL_GENOMES_MODULES + "/" + str(ncbi_id_donor) + "/*.json"):

      with open(gfile, 'r') as f:
         donors_genome_mos = json.load(f)

      counter_beneficiary = 0
      for missing_parts in list_of_missing_parts:

         # BENEFICIARY GENOME: beneficiayrs_genome_files[counter_beneficiary]
         genome_pair = "_".join([str(counter_donor), str(counter_beneficiary)])

         pair_metabolic_interactions[genome_pair] = {}
         pair_metabolic_interactions[genome_pair]["beneficiary_genome"] = beneficiayrs_genome_files[counter_beneficiary]
         pair_metabolic_interactions[genome_pair]["doner_genome"]       = gfile.split("/")[-1]
         pair_metabolic_interactions[genome_pair]['potent-met-inter']   = {}

         for module, alternative_paths in missing_parts.items():
            if module not in donors_genome_mos.keys():
               continue

            donors_module = list(donors_genome_mos[module].values())
            donors_module = [x[3:] for x in donors_module]

            for path, missing_kos in alternative_paths.items():
               check =  all(item in donors_module for item in missing_kos)

               if check:
                  # Index here denotes the several ways that complementarity might occur in a module using 2 specific genomes 
                  if module not in pair_metabolic_interactions[genome_pair]['potent-met-inter']:
                     pair_metabolic_interactions[genome_pair]['potent-met-inter'][module] = {}

                  index = len(pair_metabolic_interactions[genome_pair]['potent-met-inter'][module])
                  pair_metabolic_interactions[genome_pair]['potent-met-inter'][module][index] = {}
                  pair_metabolic_interactions[genome_pair]['potent-met-inter'][module][index]["path"] = path[1:-1].replace("'","").replace(" ", "").split(",")
                  pair_metabolic_interactions[genome_pair]['potent-met-inter'][module][index]["complements"] = missing_kos

         counter_beneficiary += 1
      counter_donor += 1

   return pair_metabolic_interactions

def set_default(obj):
    if isinstance(obj, set):
        return list(obj)
    raise TypeError

def pathway_complementarity(beneficiary_ncbi_id, doner_ncbi_id):

   # Load KEGG modules definitions
   definitions       = open(KMODULES_DEFINITIONS_PARSED, "r")
   definitions_map   = json.load(definitions)

   gaps, gfiles      = extract_alternatives_to_gapfill(beneficiary_ncbi_id, definitions_map)
   complementarities = check_for_complements(doner_ncbi_id, gaps, beneficiary_ncbi_id, gfiles)
   pair_ids = "_".join([str(beneficiary_ncbi_id), str(doner_ncbi_id)])
   namefile = ".".join([pair_ids, "json"])
   out_file = open(namefile, "w")
   json.dump(complementarities, out_file, default=set_default)

   return complementarities

if __name__ == "__main__": 
   """
   Run an example case for how this module works.
   Use KEGG genomes so you can easily display those as KEGG maps 
   """
   pathway_complementarity(1520, 476)
   pathway_complementarity(476, 1520)
