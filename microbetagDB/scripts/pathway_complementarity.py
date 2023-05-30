#!/usr/bin/env python
"""
license      = "GPLv3"
description  = 
author       = Haris Zafeiropoulos
author_email = haris.zafeiropoulos@kuleuven.be
"""
import sys, json, os
import glob
import itertools

microbetagDB_path                = "/".join( os.path.abspath(__file__).split("/")[:-2] )
KEGG_MODULES_DEFINITIONS_PARSED  = os.path.join(microbetagDB_path, "mappings/kegg_mappings/module_definition_map.json")
ALL_GENOMES_MODULES              = os.path.join(microbetagDB_path,"ref-dbs/all_genome_modules")
structurals = [ "md:M00144","md:M00149","md:M00151","md:M00152","md:M00154","md:M00155","md:M00153", "md:M00156", "md:M00158", "md:M00160" ]

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

def set_default(obj):
   """
   Converts a set to a list. If obj not a set, then raises an error. 
   """
   if isinstance(obj, set):
      return list(obj)
   raise TypeError

def build_genome_url(genomeId, genome_filename):
   """
   Returns the link to the genome under study
   """
   if genomeId[:2] == "GC":
      genomeId = "_".join(genome_filename.split("_")[:2])
      url = "https://www.ncbi.nlm.nih.gov/datasets/genome/" + genomeId
   elif genomeId[:4] == "MGYG":
      url = "https://www.ebi.ac.uk/metagenomics/genomes/" + genomeId
   else:
      url = "https://www.kegg.jp/kegg-bin/show_organism?org=" + genomeId

   return url, genomeId


def extract_alternatives_to_gapfill(ncbi_id_beneficiary, mo_map):
   """
   This function aims at finding KO terms missing for the beneficiary species, that if supplied, 
   could complete a module, based on all the possible alternative paths based on the module's definition. 

   definition_under_study_proc: a list of lists, where any inner list needs to be fully covered in order the step to be considered as legit. 
                                e.g. [['K00134', 'K00927'], ['K00150', 'K00927'], ['K11389']] 
                                 means thath we need either both ['K00134', 'K00927'] or both ['K00150', 'K00927'] or just 'K11389'.

   """
   list_of_alternatives_to_gap = []
   beneficiarys_genomes = []

   # [TODO] Consider how we will handle the unity of the alternatives_to_gap
   for gfile in glob.glob(ALL_GENOMES_MODULES + "/" + str(ncbi_id_beneficiary) + "/*.json"):

      # Load a genome of the beneficiary NCBI TaxID with the KO terms found on it related to each module
      with open(gfile, 'r') as f:
         genome_mos = json.load(f)

      # Species in need of by-products
      complete_modules    = []
      alternatives_to_gap = {}

      # Parse genome's KEGG modules and terms as found in the KEGG MODULES db
      counter_complete_modules = 0
      for module, kos_on_its_own in genome_mos.items(): 

         # Skip structural modules
         if module in structurals:
            continue

         # Get KOs related to the module under study that are present on the beneficiary's genome
         list_of_kos_present = [ k for k in kos_on_its_own.values() ]

         # Get the definition of the module under study 
         definition_under_study      = mo_map[module]['steps']
         definition_under_study_proc = []
         for term in definition_under_study.values():
            if isinstance(term, str):
               term = [term]
               print("HOW IS THIS POSSIBLE ? ")
            definition_under_study_proc.append(term)

         # Get all the possible combinations to have a complete path 
         potential_compl_paths =  [list(tup) for tup in itertools.product(*definition_under_study_proc)]

         # Edit these paths so they are in a flattened format
         flat_potent_compl_paths = []

         for path in potential_compl_paths:
            flat_potent_compl_paths.append(flatten(path))

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
                  # alternatives_to_gap[module]["1"] = {}                  
               else:
                  alternatives_to_gap[module][str(path)] = gaps
               #    alternatives_to_gap[module][str( len(alternatives_to_gap[module]) + 1 )] = {}
               # alternatives_to_gap[module][str(len(alternatives_to_gap[module]))]["path"] = str(path)
               # alternatives_to_gap[module][str(len(alternatives_to_gap[module]))]["gaps"] = gaps

      # Remove any complete module
      for key in complete_modules:
         if key in alternatives_to_gap:
            del alternatives_to_gap[key]


      # Keep only the most feasible gap fills
      # [ATTENTION] This part is the one requiring for the most computational time
      for module, path_gaps in alternatives_to_gap.items():

         tmp  = tmp2 = alternatives_to_gap[module].copy()
         min_val = min([len(path_gaps[ele]) for ele in path_gaps])

         # Make sure we do not keep cases that are supersets of another
         shortest_alternatives = [  list(tmp2.keys())[ list(tmp2.values()).index(s) ]  for s in tmp2.values() if not any( s.issuperset(i) and len(s) > len(i) for i in tmp2.values() ) ]
         """
         [ATTENTION] Here we remove potential cases of complementarity IF the KO terms required are more than 1 than the minimum required OR if they are not included in the shortest_alternatives
         """
         for path, gaps in alternatives_to_gap[module].items():
            if len(gaps) > min_val + 1 or path not in shortest_alternatives:
               del tmp[path]

         alternatives_to_gap[module] = tmp

      list_of_alternatives_to_gap.append(alternatives_to_gap)
      beneficiarys_genomes.append(gfile.split("/")[-1])

      print("NCBI Tax id:", ncbi_id_beneficiary, "\tbeneficiary genome:", gfile.split("/")[-1], "\t# of complete modules:\t", str(counter_complete_modules))
      print("~~~")


   return list_of_alternatives_to_gap, beneficiarys_genomes

def check_for_complements(ncbi_id_donor, list_of_complPath_gaps, beneficiary_ncbi_id, beneficiayrs_genome_files):
   """ 
   Check whether the potential giver (benefactor), can gap fill one of the alternatives 
   extracted by the extract_alternatives_to_gapfill() function
   list_of_complPath_gaps:  a list where each entry corresponds to a beneficiary genome and consists of a dictionary 
                           where keys are the KEGG MODULE ids and values, a dictionary with a complete alternative path of the module as key
                           and value its missing KO terms. E.g.   'md:M00881': {"['K03801', 'K03644']": {'K03801'}}
   """
   counter_donor = 0
   pair_metabolic_interactions = {}
   # pair_metabolic_interactions['beneficiary'] = {}
   # pair_metabolic_interactions['donor'] = {}
   # pair_metabolic_interactions['beneficiary']['species_name'] = 
   # pair_metabolic_interactions['beneficiary']['NCBI_ID']      = "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=" + str(beneficiary_ncbi_id)
   # pair_metabolic_interactions['donor']['species_name']       = 
   # pair_metabolic_interactions['donor']['NCBI_ID']            ="https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=" + str(ncbi_id_donor)

   pair_metabolic_interactions['beneficiary-NCBI-Id'] = str(beneficiary_ncbi_id)
   pair_metabolic_interactions['donor-NCBI-Id']       = str(ncbi_id_donor)

   for gfile in glob.glob(ALL_GENOMES_MODULES + "/" + str(ncbi_id_donor) + "/*.json"):

      with open(gfile, 'r') as f:
         donors_genome_mos = json.load(f)

      counter_beneficiary = 0
      for complPath_gaps in list_of_complPath_gaps:

         # BENEFICIARY GENOME: beneficiayrs_genome_files[counter_beneficiary]
         genome_pair = "_".join([str(counter_donor), str(counter_beneficiary)])

         # ------------------------------------
         # Get urls of the genomes involved
         beneficiary_genome_filename           = beneficiayrs_genome_files[counter_beneficiary]
         beneficiary_genomeId                  = beneficiary_genome_filename.split("_")[0]
         beneficiary_url, beneficiary_genomeId = build_genome_url(beneficiary_genomeId, beneficiary_genome_filename)

         donor_genome_filename       = gfile.split("/")[-1]
         donors_genomeId             = donor_genome_filename.split("_")[0]
         donors_url, donors_genomeId = build_genome_url(donors_genomeId, donor_genome_filename)
         # -----------------------------------

         pair_metabolic_interactions[genome_pair] = {}
         pair_metabolic_interactions[genome_pair]["beneficiary_genome"] = beneficiary_url


         pair_metabolic_interactions[genome_pair]["donor_genome"]       = donors_url
         pair_metabolic_interactions[genome_pair]['potent-met-inter']   = {}


         for module, alternative_paths in complPath_gaps.items():

            # In case a module is not present at all in the donor, i.e. the donor has no KO that's part of this module, go to the next one
            if module not in donors_genome_mos.keys():
               continue

            # Keep a list with all the KOs on the donor related to the module under study
            donors_module = list(donors_genome_mos[module].values())

            complements = []
            for path, missing_kos in alternative_paths.items():
               
               # Check whether the donor has all the KOs necessary according to a specific path
               check =  all(item in donors_module for item in missing_kos)

               if check:

                  # Index here denotes the several ways that complementarity might occur with different KOs provided 
                  if module not in pair_metabolic_interactions[genome_pair]['potent-met-inter']:
                     pair_metabolic_interactions[genome_pair]['potent-met-inter'][module] = {}

                  if missing_kos not in complements:
                     pc_index = len(pair_metabolic_interactions[genome_pair]['potent-met-inter'][module])
                     pair_metabolic_interactions[genome_pair]['potent-met-inter'][module][pc_index] = {}
                     pair_metabolic_interactions[genome_pair]['potent-met-inter'][module][pc_index]["complement"] = missing_kos
                     pair_metabolic_interactions[genome_pair]['potent-met-inter'][module][pc_index]["path"] = path[1:-1].replace("'","").replace(" ", "").split(",")
                     complements.append(missing_kos)
                  else:
                     pc_index   = complements.index(missing_kos)
                     path_index = len(pair_metabolic_interactions[genome_pair]['potent-met-inter'][module][pc_index]["path"])
                     pair_metabolic_interactions[genome_pair]['potent-met-inter'][module][pc_index]["path"][path_index] = path[1:-1].replace("'","").replace(" ", "").split(",")                  

         counter_beneficiary += 1
      counter_donor += 1

   return pair_metabolic_interactions

def pathway_complementarity(beneficiary_ncbi_id, donor_ncbi_id):
   print("\n\n ~~~ PATHWAY COMPLEMENTARITY PRECALCULATION STEP ~~~")
   # Load KEGG modules definitions
   definitions       = open(KEGG_MODULES_DEFINITIONS_PARSED, "r")
   definitions_map   = json.load(definitions)

   print("Extract alternatives that could be gapfilled")
   gaps, gfiles      = extract_alternatives_to_gapfill(beneficiary_ncbi_id, definitions_map)

   print("\n\n *** \nGet potential complements from the donor")
   complementarities = check_for_complements(donor_ncbi_id, gaps, beneficiary_ncbi_id, gfiles)

   pair_ids = "_".join([str(beneficiary_ncbi_id), str(donor_ncbi_id)])
   namefile = ".".join([pair_ids, "json"])
   out_file = open(namefile, "w")

   json.dump(complementarities, out_file, default = set_default)
   return complementarities

if __name__ == "__main__": 
   """
   Run an example case for how this module works.
   Use KEGG genomes so you can easily display those as KEGG maps 
   """
   pathway_complementarity(476, 1520)

