"""
Parse user's input OTU table
"""
import pandas as pd
import logging
import os, sys
from utils.variables import * 

def count_comment_lines(my_otu_table, my_taxonomy_column):

   skip_rows = 0
   with open(my_otu_table, 'r') as f:
      for line in f:
         if line.startswith('#') and my_taxonomy_column not in line :
            skip_rows += 1
         elif my_taxonomy_column in line: 
            line     = line.split("\t")
            line[-1] = line[-1][:-1]
         else:
               break
   return skip_rows


def is_tab_separated(my_otu_table, my_taxonomy_column, otu_identifier_column):

   number_of_commented_lines = count_comment_lines(my_otu_table, my_taxonomy_column)

   otu_table = pd.read_csv(my_otu_table, sep = "\t", skiprows= number_of_commented_lines)

   if otu_table.shape[1] < 2:

      logging.error("The OTU table provided is not a tab separated file. Please convert your OTU table to .tsv or .csv format.")

   else: 

      try: 

         # Build a 2-col dataframe (identifier - taxonomy)
         taxonomies = otu_table.filter([otu_identifier_column, my_taxonomy_column])

         # Split the taxonomy column and split it based on semicolumn! -----  MAKE THAT MORE GENERAL AT A LATER POINT
         splitted_taxonomies         = taxonomies[my_taxonomy_column].str.split(';', expand = True)
         splitted_taxonomies.columns = ["Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
         splitted_taxonomies[otu_identifier_column] = taxonomies[otu_identifier_column]

         # Make sure there is no white space character in the Species, Genus and Family levels
         splitted_taxonomies['Species'] = splitted_taxonomies['Species'].str.strip()
         splitted_taxonomies['Genus']   = splitted_taxonomies['Genus'].str.strip()
         splitted_taxonomies['Family']  = splitted_taxonomies['Family'].str.strip()

      except:
         logging.error("The taxonomy column could not be retrieved from you OTU table. \n \
                        Check for any strange characters on your OTU table or its format.")

      # Read the species, genus and family files of Silva db with their NCBI IDs
      silva_species_ncbi_id         = pd.read_csv(SILVA_SPECIES_NCBI_ID, sep = "\t")
      silva_species_ncbi_id.columns = ['Species', 'ncbi_tax_id']
      silva_species_ncbi_id['Species'].str.strip()

      silva_genus_ncbi_id         = pd.read_csv(SILVA_GENUS_NCBI_ID, sep = "\t")
      silva_genus_ncbi_id.columns = ['Genus', 'ncbi_tax_id']
      silva_genus_ncbi_id['Genus'].str.strip()

      silva_family_ncbi_id         = pd.read_csv(SILVA_FAMILY_NCBI_ID, sep = "\t")
      silva_family_ncbi_id.columns = ['Family', 'ncbi_tax_id']
      silva_family_ncbi_id['Family'].str.strip()

      # Build a dataframe for the Species, Genus and Family taxonomies present on the OTU table along with their corresponding NCBI Tax IDs
      species_present  = silva_species_ncbi_id.merge(splitted_taxonomies, on = ['Species'])
      species_present  = species_present[['Species', 'ncbi_tax_id', otu_identifier_column]]
      genus_present    = silva_genus_ncbi_id.merge(splitted_taxonomies, on = ['Genus'])
      families_present = silva_family_ncbi_id.merge(splitted_taxonomies, on = ['Family'])


   return otu_table, species_present, genus_present, families_present


def ensure_flashweave_format(my_otu_table, my_taxonomy_column, otu_identifier_column):

   flashweave_table = my_otu_table.drop(my_taxonomy_column, axis = 1)
   float_col        = flashweave_table.select_dtypes(include=['float64']) 
   
   for col in float_col.columns.values:
      flashweave_table[col] = flashweave_table[col].astype('int64')

   my_otu_table[otu_identifier_column] = flashweave_table[otu_identifier_column]
   file_to_save                  = os.path.join(FLASHWEAVE_OUTPUT_DIR, "otu_table_flashweave_format.tsv")

   flashweave_table.to_csv(file_to_save, sep ='\t', index = False)

   return my_otu_table


def edgelist_to_ncbi_ids(edgelist, species_to_ncbi_ids, otu_identifier_column):


   if EDGE_LIST: 
      associations = pd.read_csv(edgelist, sep = "\t", skiprows = 2)

   else: 

      associations = pd.read_csv(edgelist, sep = "\t", skiprows = 2)


   associations.columns = ['taxon_A', 'taxon_B', 'evidence']


   left_part = associations[['taxon_A', 'evidence']]
   left_part.columns = [otu_identifier_column, 'evidence']
   associations_with_ncbiIds_left = left_part.merge(species_to_ncbi_ids, on = otu_identifier_column)
   associations_with_ncbiIds_left.columns = [otu_identifier_column, 'evidence', 'Taxon_A', 'NCBI_ID_A' ]


   right_part = associations[['taxon_B', 'evidence']]
   right_part.columns = [otu_identifier_column, 'evidence']
   associations_with_ncbiIds_right = right_part.merge(species_to_ncbi_ids, on = otu_identifier_column)
   associations_with_ncbiIds_right.columns = [otu_identifier_column, 'evidence', 'Taxon_B', 'NCBI_ID_B' ]
   associations_with_ncbiIds_right.drop("evidence", axis=1, inplace=True)


   data = [associations_with_ncbiIds_left['Taxon_A'], associations_with_ncbiIds_left['NCBI_ID_A'], associations_with_ncbiIds_right['Taxon_B'], associations_with_ncbiIds_right['NCBI_ID_B>

   headers = ["taxon_a", "ncbi_a", "taxon_b", "ncbi_b", "evidence"]
   
   associations_with_ncbiIds = pd.concat(data, axis = 1, keys = headers)
   
   print(associations_with_ncbiIds)
   print("\n\n~~~~~~\n\n")
   print(associations)   

   return associated_pairs



