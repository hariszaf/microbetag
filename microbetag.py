#!/usr/bin/env python 

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
   Assure the output directory
   """
   if not os.path.exists(OUT_DIR):
      os.mkdir(OUT_DIR)

   """
   Make sure OTU table in tab separated format
   """
   print(TAX_COL)
   is_tab_separated(OTU_TABLE, TAX_COL)


   """
   STEP: FlashWeave 
   """
   if not args.el:

      logging.info("Assure OTU table format fits FlashWeave")
      if not os.path.exists(FLASHWEAVE_OUTPUT_DIR):
         os.mkdir(FLASHWEAVE_OUTPUT_DIR)
      ensure_flashweave_format(OTU_TABLE, TAX_COL, OTU_COL)
      
      logging.info("Run FlashWeave")


      flashweave_parmas = [
         "julia", FLASHWEAVE_SCRIPT, FLASHWEAVE_TMP_INPUT
      ]

      flashweave_command = ' '.join(flashweave_parmas)
      os.system(flashweave_command)
      
      sys.exit(0)




   """
   STEP: FAPROTAX
   """
   logging.info('STEP: FAPROTAX database oriented analaysis'.center(50, '*'))
   if args.i: 

      faprotax_check = False

      # if args.t is None:
      #    with open(OTU_TABLE) as f:
      #       lines = f.readlines()
      #    lengths = set()
      #    for line in lines: 
      #       line = line.split("\t")
      #       lengths.add(len(line))
      #    max_num_of_cells = max(lengths)
      #    for line in lines: 
      #       line = line.split("\t")
      #       if len(line) == max_num_of_cells:
      #          TAX_COL = '"' + line[-1].rstrip() + '"'
      #          OTU_COL = line[0].rstrip()
      #          break
      #       if args.com is None:
      #          COM_CHAR = '"' + line[0][0] + '"'

      faprotax_params = [
         "python", FAPROTAX_SCRIPT,
         "-i", OTU_TABLE,
         "-o", FAPROTAX_FUNCT_TABLE,
         "-g", FAPROTAX_DB,
         "-c", COM_CHAR,
         "-d", '"' + TAX_COL + '"',    # "-d", TAX_COL,
         "-v",
         "-s", FAPROTAX_SUB_TABLES,
         "-force"
      ]

      if args.ch:
         faprotax_params = faprotax_params + ["--column_names_are_in", COM_HEAD]

      cmd = ' '.join(faprotax_params)

      try:
         logging.info('Phenotypic analysis using BugBase')
         logging.info(cmd)
         os.system(cmd)
         faprotax_check = True
      except:
         logging.exception("\nSomething went wrong when running the BugBase analysis!")

      if faprotax_check: 
         path_to_subtables = os.path.join(BASE, FAPROTAX_SUB_TABLES)
         print(">>>>> PATH; ", path_to_subtables)
         otu_faprotax_functions_assignment(path_to_subtables)






if __name__ == '__main__':
    main()