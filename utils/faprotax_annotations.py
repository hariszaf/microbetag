"""
Parse the sub tables of the faprotax analysis 
to assign the biological processes related to each OTU 
"""

__author__  = 'Haris Zafeiropoulos'
__email__   = 'haris-zaf@hcmr.gr'
__status__  = 'Development'
__license__ = 'GPLv3'

import os

def otu_faprotax_functions_assignment(path_to_subtables):

   otu_faprotax_assignments = {}

   for process_name in os.listdir(path_to_subtables): 
      
      f = os.path.join(path_to_subtables, process_name)
      table_file = open(f, "r")
      table_file = table_file.readlines()

      for line in table_file[2:]:

         otu_id = line.split("\t")[1]

         if otu_id not in otu_faprotax_assignments:

            otu_faprotax_assignments[otu_id] = [process_name]
         
         else:

            otu_faprotax_assignments[otu_id].append(process_name)

   return otu_faprotax_assignments


         


   






