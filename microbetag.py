#!/usr/bin/env python3

__author__  = 'Haris Zafeiropoulos'
__email__   = 'haris-zaf@hcmr.gr'
__status__  = 'Development'
__license__ = 'GPLv3'


from utils import *
import os, logging


def main():

   """
   Setting logging
   """
   # Set logger
   logger = logging.getLogger()
   logger.setLevel(logging.INFO)
   formatter = logging.Formatter(
      '%(asctime)s [%(levelname)-8s] %(message)s',
      datefmt='%Y-%m-%d %H:%M:%S')


   """
   Assure the output directory
   """
   if not os.path.exists(OUT_DIR):
      os.mkdir(OUT_DIR)


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
   logging.info('User input: {}'.format(' '.join(sys.argv)))





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
      species_present     = get_species(ext, TAX_COL, OTU_COL)
      # # Get a list of dictionaries with the OTU ID, the species name, the microbetag id and the NCBI Taxonomy id for each entry 
      # # e.g.
      # # {'#OTU ID': 821, 'taxonomy': 'D_6__planctomycete str. 394', 'microbetag_id': 'microbetag_821', 'ncbi_tax_id': '79876'}
      # otu_species_ncbi_id = species_present.to_dict(orient = 'records')

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
      species_pairs = edge_list_of_ncbi_ids(FLASHWEAVE_EDGELIST, species_present)



   print(species_pairs)
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