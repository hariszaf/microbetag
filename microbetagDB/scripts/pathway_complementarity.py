"""
license      = "GPLv3"
description  = ""
author       = Haris Zafeiropoulos
author_email = haris.zafeiropoulos@kuleuven.be
"""
import sys, json, os
import glob, csv
import itertools
import multiprocessing
import time

microbetagDB_path                = "/".join( os.path.abspath(__file__).split("/")[:-2] )
KEGG_MODULES_DEFINITIONS_PARSED  = os.path.join(microbetagDB_path, "mappings/kegg_mappings/module_definition_map.json")
ALL_GENOMES_MODULES_PATH         = os.path.join(microbetagDB_path, "ref-dbs/all_genome_modules")
COMPLEMENTARITIES_PATH           = os.path.join(microbetagDB_path, "ref-dbs/complementarities")
KO_ANNOTATED_GENOMES_FILE        = os.path.join(microbetagDB_path, "ref-dbs/ko_per_mod_in_genome_per_ncbiId.json")
structurals = [ "md:M00144","md:M00149","md:M00151","md:M00152","md:M00154","md:M00155","md:M00153", "md:M00156", "md:M00158", "md:M00160" ]


def pathway_complementarity(beneficiary_ncbi_id):
   """
   This function is going to parse all the other ncbi id related genomes to check for complements
   for the species under study
   """

   # Check if beneficiary ncbi id has already been parsed
   if beneficiary_ncbi_id in ncbiIds_studied:
      processed_ncbi_ids.remove(beneficiary_ncbi_id)
      print("NCBI id already parsed my friend.")
      return processed_ncbi_ids

   # Get all the possible alternatives for the KEGG modules of the beneficiary to gapfilled
   gaps, gfiles      = extract_alternatives_to_gapfill(beneficiary_ncbi_id, definitions_map)

   if gaps == 0:
      return 1

   logfile  = open("logfile", "a")

   # Get potential gapfills from all the rest ncbiIds
   donors_ncbiIds = ncbi_ids.copy()
   donors_ncbiIds.remove(beneficiary_ncbi_id)
   for donor_ncbi_id in donors_ncbiIds:
      # Get complmementarities for specific pair
      check_for_complements(donor_ncbi_id, gaps, beneficiary_ncbi_id, gfiles)

   # Keep track of the pairs explored
   logfile.write(beneficiary_ncbi_id + "\n")

   # Remove the pair from the running species
   try:
      processed_ncbi_ids.remove(beneficiary_ncbi_id)
   except:
      print("you must be in your last chunk. otherwise you have an error... doom doom dooom")

   return processed_ncbi_ids


def chunker(seq, size):
   """
   get chunks from a list based on a size
   """
   return (seq[pos:pos + size] for pos in range(0, len(seq), size))


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


def build_kegg_url():
   """
   Build url to colorify the related to the module kegg map based on the KO terms
   of the beneficiary (pink) and those it gets from the donor (green)
   """
   # Load the dictionary with the kegg modules and their corresponding maps
   color_mapp_base_url = "https://www.kegg.jp/kegg-bin/show_pathway?"
   present_kos_color   = "%09%23EAD1DC/"
   complemet_kos_color = "%09%2300A898/"
   maps = open(
            os.path.join(microbetagDB_path, "mappings/kegg_mappings/module_map_pairs.tsv"),
            "r"
         )
   module_map = {}
   for line in maps:
      module, mmap = line.split("\t")
      module_map[module[:-1]] = mmap[1:-1]

   # Make a url pointing at a colored kegg map based on what's on the beneficiary's genome and what it gets as complement from the donor
   beneficiarys_kos = ""
   complements_kos  = ""
   for ko_term in clean_path:
      if ko_term not in missing_kos:
         beneficiarys_kos = "".join([beneficiarys_kos, ko_term, present_kos_color]) # += ko_term + present_kos_color
      else:
         complements_kos = "".join([beneficiarys_kos, ko_term, complemet_kos_color])  # += ko_term + complemet_kos_color

   try:
      """ we do the try as in some rare cases, the module might not have a related map"""
      url_ko_map_colored = "".join([color_mapp_base_url, module_map[module],  "/", beneficiarys_kos, complements_kos])
   except:
      url_ko_map_colored = "N/A"

   return url_ko_map_colored


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

   try:
      test = ko_annotated_genomes_per_ncbiId[ncbi_id_beneficiary]
   except:
      print( "No genomes for beneficiary: ",  ncbi_id_beneficiary,)
      return 0, 0

   for beneficary_genome_file, beneficiarys_annotated_genomes in ko_annotated_genomes_per_ncbiId[ncbi_id_beneficiary].items():

      # Species in need of by-products
      complete_modules    = set()
      alternatives_to_gap = {}

      # Parse genome's KEGG modules and terms as found in the KEGG MODULES db
      counter_complete_modules = 0
      for module, kos_on_its_own in beneficiarys_annotated_genomes.items():

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
               sys.exit(0)
            definition_under_study_proc.append(term)

         # Get all the possible combinations to have a complete path
         potential_compl_paths =  [ list(tup) for tup in itertools.product(*definition_under_study_proc) ]

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
                  complete_modules.add(module)
            else:
               """Keep track of the missing steps in the alternative
               alteritves_to_gap will look like this:
               'md:M00550': {"['K02821', 'K02822', 'K03077', 'K03078', 'K03079', 'K03475', 'K03476']": {'K03079', 'K03476', 'K03078'}}
               where the values of a module is a dictionary
               where a path is the key and the missing KOs from the genome is the value"""

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
      # [ATTENTION] This part is the one requiring for the most computational time
      for module, path_gaps in alternatives_to_gap.items():

         tmp  = tmp2 = alternatives_to_gap[module].copy()
         min_val = min([len(path_gaps[ele]) for ele in path_gaps])

         # Make sure we do not keep cases that are supersets of another
         """ we make this instead of the muted in terms of faster implementation.. """
         values = list(tmp2.values())
         shortest_alternatives = [list(tmp2.keys())[values.index(s)] for s in values if not any(s.issuperset(i) and len(s) > len(i) for i in values)]
         # shortest_alternatives = [  list( tmp2.keys() )[ list(tmp2.values()).index(s) ]  for s in tmp2.values() if not any( s.issuperset(i) and len(s) > len(i) for i in tmp2.values() ) ]

         for path, gaps in alternatives_to_gap[module].items():
            """
            [ATTENTION] Here we remove potential cases of complementarity IF the KO terms required are more than 1 than the minimum required OR if they are not included in the shortest_alternatives
            """
            if len(gaps) > min_val + 1 or path not in shortest_alternatives:
               del tmp[path]

         alternatives_to_gap[module] = tmp

      list_of_alternatives_to_gap.append(alternatives_to_gap)
      beneficiarys_genomes.append(beneficary_genome_file)

   return list_of_alternatives_to_gap, beneficiarys_genomes


def check_for_complements(donor_ncbi_id, set_of_complPath_gaps, beneficiary_ncbi_id, beneficiayrs_genome_files):
   """
   Check whether the potential giver (benefactor), can gap fill one of the alternatives
   extracted by the extract_alternatives_to_gapfill() function
   set_of_complPath_gaps:  a set where each entry corresponds to a beneficiary genome and consists of a dictionary
                           where keys are the KEGG MODULE ids and values, a dictionary with a complete alternative path of the module as key
                           and value its missing KO terms. E.g.   'md:M00881': {"['K03801', 'K03644']": {'K03801'}}
   """
   counter_donor = 0
   tab_file = open("/staging/leuven/stg_00106/haris/microbetag/COMPLEMENTARITIES/complementarities_" + str(beneficiary_ncbi_id)  + ".tsv", "a+")
   # Parse corresponding genomes of the donor's NCBI Tax Id
   for donors_genome_file in glob.glob(ALL_GENOMES_MODULES_PATH + "/" + str(donor_ncbi_id) + "/*.json"):

      with open(donors_genome_file, 'r') as f:
         donors_genome_mos = json.load(f)

      counter_beneficiary = 0
      for complPath_gaps in set_of_complPath_gaps:

         # Get urls of the genomes involved
         beneficiary_genome_filename           = beneficiayrs_genome_files[counter_beneficiary]
         beneficiary_genomeId                  = beneficiary_genome_filename.split("/")[-1].replace("_kos_related_to_mos.json", "")

         #donor_genome_filename       = donors_genome_file.split("/")[-1]
         donors_genomeId             = donors_genome_file.split("/")[-1].replace("_kos_related_to_mos.json", "")

         for module, alternative_paths in complPath_gaps.items():

            # In case a module is not present at all in the donor, i.e. the donor has no KO that's part of this module, go to the next one
            if module not in donors_genome_mos.keys():
               continue

            # Keep a list with all the KOs on the donor related to the module under study
            donors_module = list(donors_genome_mos[module].values())
            for path, missing_kos in alternative_paths.items():

               clean_path = path[1:-1].replace("'","").replace(" ", "").split(",")

               # Check whether the donor has all the KOs necessary according to a specific path
               check =  all(item in donors_module for item in missing_kos)

               if check:

                  # Add entry in the tab file - each entry corresponds to a complement
                  #with open("complementarities.tsv", "a") as tab_file:
                  tab_file_writer = csv.writer(tab_file, delimiter = "\t",  quoting = csv.QUOTE_MINIMAL)
                  tab_file_writer.writerow([
                     beneficiary_genomeId,
                     donors_genomeId,
                     beneficiary_ncbi_id,
                     donor_ncbi_id, module,
                     ";".join(missing_kos),
                     ";".join(clean_path)
                  ])


         counter_beneficiary += 1
      counter_donor += 1

   return 1



if __name__ == "__main__":
   """
   Run an example case for how this module works.
   Use KEGG genomes so you can easily display those as KEGG maps
   """

   # Open files to append to
   logfile  = open("logfile", "r")
   #tab_file = open("complementarities.tsv", "a")

   # Load KEGG modules definitions
   definitions       = open(KEGG_MODULES_DEFINITIONS_PARSED, "r")
   definitions_map   = json.load(definitions)

   # Load all annotated genomes per ncbiId
   ko_annotated_genomes_per_ncbiId = open(KO_ANNOTATED_GENOMES_FILE, "r")
   ko_annotated_genomes_per_ncbiId = json.load(ko_annotated_genomes_per_ncbiId)

   # Get all possible NCBI Id pairs
   ncbi_ids = set()
   for i in os.listdir(ALL_GENOMES_MODULES_PATH):
      ncbi_ids.add(str(i))

   # Get pairs already parsed in previous runs; ATTENTION! if you have new genomes for a ncbiId this will skip those
   logfile = logfile.readlines()
   ncbiIds_studied = set()
   for line in logfile:
      ncbiIds_studied.add(line[:-1])

   ncbi_ids = ncbi_ids - ncbiIds_studied

   # Set the maximum number of processes that our script is going to use
   number_processes = int(sys.argv[1])
   ncbi_ids = list(ncbi_ids)
   processed_ncbi_ids = ncbi_ids.copy()
   ncbi_id_chunks = list(chunker(ncbi_ids, number_processes))

   while len(processed_ncbi_ids) > 0:

      counter = 1
      # for index, a in enumerate(ncbi_ids):
      for chunk in ncbi_id_chunks:

         remaining_ncbiIds = open("remainingIds.tsv", "w")
         if isinstance(processed_ncbi_ids[0], list):
            processed_ncbi_ids = list(set(processed_ncbi_ids[0]).intersection(*processed_ncbi_ids[1:]))
         for line in processed_ncbi_ids:
             remaining_ncbiIds.write(line + "\n")
         remaining_ncbiIds.flush()
         time.sleep(1)

         print("\n\n**** species to go:", str(len(processed_ncbi_ids)), "  *****")
         print("**** beneficiary under study:", str(chunk), "  *****\n\n")

         # Start a pool
         pool = multiprocessing.Pool(number_processes)

         # Highly important line:
         processed_ncbi_ids = pool.map(pathway_complementarity, chunk)

         pool.close()
         pool.join()


   """in case you want to try this only for a pair of ncbi ids, mute the while loop and unmute the following lines
   you can then replace the ncbi ids with those of your choice"""

   # pairs_of_ncbi_ids_for_ncbiId_a = [("476", "1520"), ("1520", "476")]
   # pathway_complementarity(("476", "1520"))
   # pathway_complementarity(("1520", "476"))


