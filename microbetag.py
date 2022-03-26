#!/usr/bin/env python3

"""
what is microbetag about and how to use it 
"""

__author__  = 'Haris Zafeiropoulos'
__email__   = 'haris-zaf@hcmr.gr'
__status__  = 'Development'
__license__ = 'GPLv3'
__version__ = 'v.0.0.1'

from utils import *
import os

def main():

   """
   Assure the output directory
   """
   if not os.path.exists(OUT_DIR):
      os.mkdir(OUT_DIR)

   """
   Setting logging
   """
   # Using FileHandler writing log to file
   logfile = os.path.join(OUT_DIR, 'log.txt')
   fh      = logging.FileHandler(logfile)
   fh.setLevel(logging.DEBUG)
   fh.setFormatter(formatter)

   # Using StreamHandler writing to console
   ch = logging.StreamHandler()
   ch.setLevel(logging.INFO)
   ch.setFormatter(formatter)

   # Add the two Handlers
   logger.addHandler(ch)
   logger.addHandler(fh)

   """
   Welcome message and arguments values
   """
   logging.info("Hello microbe-fun! microbetag is about to start!")
   logging.info('Your command was: {}'.format(' '.join(sys.argv)))


   """
   STEP: OTU table preprocess 
   """
   if OTU_TABLE: 
      logging.info("Make sure OTU table in tab separated format")

      # Load the initial OTU table as a pandas dataframe
      otu_table = is_tab_separated(OTU_TABLE, TAX_COL)
      logging.info("Your OTU table is a tab separated file that microbetag can work with.")


      if not EDGE_LIST:
         """
         Pre-process
         """
         logging.info("The user has not provided an edge list. microbetag will build one using FlashWeaeve.")
         if not os.path.exists(FLASHWEAVE_OUTPUT_DIR):
            os.mkdir(FLASHWEAVE_OUTPUT_DIR)

         logging.info("Assure OTU table format fits FlashWeave")
         ext = ensure_flashweave_format(otu_table, TAX_COL, OTU_COL)

      else:

         ext = otu_table.copy()
         ext['microbetag_id'] = otu_table[OTU_COL]

      logging.info("Get the NCBI Taxonomy id for the taxonomies at the species level")

      # The get_species() function returns a pandas dataframe
      species_present     = get_species(ext, TAX_COL, OTU_COL)  


   """
   STEP: Get co-occurrence network
   """
   logging.info('STEP: Get co-occurrence network'.center(50, '*'))
   if not EDGE_LIST:

      """
      Run FlashWeave
      """
      logging.info("Run FlashWeave")
      flashweave_parmas = [
         "julia", FLASHWEAVE_SCRIPT, FLASHWEAVE_OUTPUT_DIR, FLASHWEAVE_TMP_INPUT
      ]
      flashweave_command = ' '.join(flashweave_parmas)
      os.system(flashweave_command)

      """
      Assign NCBI Taxonomy ids to the nodes of the network
      """
      logging.info("Match a NCBI Taxonomy id to the species present on the OTU table.")

      # The edge_list_of_ncbi_ids() function returns a dictionary like this: 
      # {'0': {taxon_1: {}, taxon_2 }
      edge_list = edge_list_of_ncbi_ids(FLASHWEAVE_EDGELIST, species_present)   


   """
   STEP: PATHWAY COMPLEMENTARITY
   """
   logging.info("STEP: Pathway complementarity module: metabolic interactions ".center(50, '*'))
   if PATHWAY_COMPLEMENTARITY == True: 
      

      """
      Remember to change the path from kegg_genomes to all_genomes once the latter is ready 
      """

      set_of_ncbi_ids_with_available_genomes = set(os.listdir('ref-dbs/kegg_genomes/'))



      if not EDGE_LIST:

         for pair in edge_list.values(): 

            taxon_a = pair['taxon_1']['ncbi_tax_id']
            taxon_b = pair['taxon_2']['ncbi_tax_id']

            if taxon_a in set_of_ncbi_ids_with_available_genomes and taxon_b in set_of_ncbi_ids_with_available_genomes: 

               pathway_complementarity(beneficiary_ncbi_id, doner_ncbi_id)

      else: 

         """
         In case you are using your own edge list or if you have already run microbetag and you have already one
         """


   sys.exit(0)



   """
   STEP: FAPROTAX
   """
   logging.info('STEP: FAPROTAX database oriented analaysis'.center(50, '*'))
   if OTU_TABLE: 

      if not os.path.exists(FAPROTAX_OUTPUT_DIR):
         os.mkdir(FAPROTAX_OUTPUT_DIR)

      faprotax_check = False

      faprotax_params = [
         "python3", FAPROTAX_SCRIPT,
         "-i",     OTU_TABLE,
         "-o",     FAPROTAX_FUNCT_TABLE,
         "-g",     FAPROTAX_DB,
         "-c", '"' + COM_CHAR + '"',
         "-d", '"' +  TAX_COL + '"',
         "-v",
         "-s",     FAPROTAX_SUB_TABLES,
      ]

      if COM_HEAD:
         faprotax_params = faprotax_params + ["--column_names_are_in", COM_HEAD]

      cmd = ' '.join(faprotax_params)

      try:
         logging.info('Phenotypic analysis using BugBase')
         logging.info(cmd)
         os.system(cmd)
         faprotax_check = True
      except:
         logging.exception("\nSomething went wrong when running the BugBase analysis!")

      # If FAPROTAX was completed, make a dictionary with OTUs as keys and processes 
      # retrieved as values
      if faprotax_check: 
         path_to_subtables = os.path.join(BASE, FAPROTAX_SUB_TABLES)
         otu_faprotax_functions_assignment(path_to_subtables)


   """
   STEP: BugBase
   """
   logging.info('STEP: BugBase database oriented analaysis'.center(50, '*'))
   if OTU_TABLE: 


      bugbase_commands = [
         "Rscript", "run.bugbase.r", 
         "-i", OTU_TABLE,
         "-o", BUGBASE_OUTPUT
      ]

      if METADATA_FILE:
         bugbase_commands = bugbase_commands + ["-m", METADATA_FILE]



      for k,v in BUGBASE_OPTS.items(): 
         if v is not None:
            print(k,v)


      cmd = ' '.join(bugbase_commands)

      try:
         logging.info('Phenotypic analysis using BugBase')
         logging.info(cmd)
         os.system(cmd)

      except:
         logging.exception("\nSomething went wrong when running the BugBase analysis!")


if __name__ == '__main__':
   main()