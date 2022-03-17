#!/usr/bin/env python 

from utils import *
import os


def main():

   logging.info('User input: {}'.format(' '.join(sys.argv)))

   """
   Welcome message
   """
   print("HELLO FRIEND")


   """
   Assure the output directory
   """
   if not os.path.exists(OUT_DIR):
      os.mkdir(OUT_DIR)


   """
   STEP: FAPROTAX
   """
   logging.info('STEP: FAPROTAX \nRun FAPROTAX analysis'.center(50, '*'))
   if args.i: 

      if args.t is None:
         with open(OTU_TABLE) as f:
            lines = f.readlines()
         lengths = set()
         for line in lines: 
            line = line.split("\t")
            lengths.add(len(line))
         max_num_of_cells = max(lengths)
         for line in lines: 
            line = line.split("\t")
            if len(line) == max_num_of_cells:
               TAX_COL = '"' + line[-1].rstrip() + '"'
               OTU_COL = line[0].rstrip()
               break
            if args.com is None:
               COM_CHAR = '"' + line[0][0] + '"'

      faprotax_params = [
         "python", FAPROTAX_SCRIPT,
         "-i", OTU_TABLE,
         "-o", FAPROTAX_OUTPUT,
         "-g", FAPROTAX_DB,
         "-c", COM_CHAR,
         "-d", TAX_COL,
         "-v"
      ]

      cmd = ' '.join(faprotax_params)
      print(cmd)
      sys.exit(0)
      try:
         logging.info('Phenotypic analysis using BugBase')
         logging.info(cmd)
         os.system(cmd)
      except:
         logging.exception("\nSomething went wrong when running the BugBase analysis!")


# collapse_table.py -i tax_table.tsv -o func_table.tsv -g FAPROTAX.txt -d "taxonomy" -c "#" -v 

!!!!   MAKE THIS AN ARGUMENT -column_names_are_in last_comment_line 

# collapse_table.py -i otu_table.tsv -o functional_table.tsv -g FAPROTAX.txt -c "#" -d "taxonomy" --omit_columns 0 --column_names_are_in last_comment_line -r report.txt -n columns_after_collapsing -v




if __name__ == '__main__':
    main()