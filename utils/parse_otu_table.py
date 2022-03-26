"""
Parse user's input OTU table
"""
import pandas as pd
import logging
import os
import sys


# Line 9 is not that clear to me.. It is my belief that the "utlis" section
# is needed cause in the package framework our root is at the /microbetag level
# If we would run this on its own, then variables would be imported without the "utils"
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


def is_tab_separated(my_otu_table, my_taxonomy_column):

   number_of_commented_lines = count_comment_lines(my_otu_table, my_taxonomy_column)

   otu_table = pd.read_csv(my_otu_table, sep = "\t", skiprows= number_of_commented_lines)

   if otu_table.shape[1] < 2:

      logging.error("The OTU table provided is not a tab separated file. Please convert your OTU table to .tsv or .csv format.")

   return otu_table


def get_species(my_otu_table, my_taxonomy_column, otu_identifier_column):


   try: 

      """
      The pd.filter() function: 
      Subset rows or columns of dataframe according to labels in the specified index. 
      Note that this routine does not filter a dataframe on its contents. 
      The filter is applied to the labels of the index.
      """

      # keep only otu ids and taxonomies assigned
      otu_id_and_taxonomy = my_otu_table.filter([otu_identifier_column, my_taxonomy_column, 'microbetag_id']) 

      # logging.info("Your OTU table lead to a dataframe with columns of the following types: \n",
                     # otu_id_and_taxonomy.dtypes.value_counts())

      # keep only the last level of lineage assigned and remove white spaces if any before the species name assigned
      otu_id_and_taxonomy[my_taxonomy_column] = otu_id_and_taxonomy[my_taxonomy_column].str.split(';').str[-1]
      otu_id_and_taxonomy[my_taxonomy_column] = otu_id_and_taxonomy[my_taxonomy_column].str.strip()

   except:

      logging.error("The taxonomy column could not be retrieved from you OTU table. \n \
                     Check for any strange characters on your OTU table or its format.")

   # Load the Silva species names along with ther NCBI Taxonomy ids
   silva_species_ncbi_id         = pd.read_csv(SILVA_SPECIES_NCBI_ID, sep = "\t")
   silva_species_ncbi_id.columns = ['species_name', 'ncbi_tax_id']
   silva_species_ncbi_id['species_name'].str.strip()
  
   # Make a column showing whether or not the species of each row is included on the Silva ref file
   otu_id_and_taxonomy['present'] = otu_id_and_taxonomy[my_taxonomy_column].isin(silva_species_ncbi_id['species_name'])

   # We make a map (dictionary) keeping as a key the column that links the 2 datadrames; in this case the species name
   map_dict = dict(zip(silva_species_ncbi_id['species_name'], silva_species_ncbi_id['ncbi_tax_id']))
   otu_id_and_taxonomy['ncbi_tax_id'] = otu_id_and_taxonomy[my_taxonomy_column].map(map_dict)

   # Keep only those rows that have "True" as a value on the 'present' column
   otu_id_species_name_ncbi_id = otu_id_and_taxonomy[otu_id_and_taxonomy.present]

   # And then remove the column 'present'
   """
   REMEMBER TO REPLACE THE FOLLOWING COMMANDS WITH .loc
   """
   otu_id_species_name_ncbi_id['ncbi_tax_id'] = otu_id_species_name_ncbi_id['ncbi_tax_id'].apply(lambda f: format(f, '.0f'))
   otu_id_species_name_ncbi_id['ncbi_tax_id'] = otu_id_species_name_ncbi_id['ncbi_tax_id'].map(str)


   logging.info("\n A table including only the taxonomies assigned to a valid species name has been built. \n \
                    The species have been mathed to their corresponding NCBI Taxonomy ids.\n")

   return otu_id_species_name_ncbi_id


def ensure_flashweave_format(my_otu_table, my_taxonomy_column, otu_identifier_column):


   flashweave_table = my_otu_table.drop(my_taxonomy_column, axis = 1)
   float_col        = flashweave_table.select_dtypes(include=['float64']) 
   
   for col in float_col.columns.values:
      flashweave_table[col] = flashweave_table[col].astype('int64')

   flashweave_table[otu_identifier_column] = 'microbetag_' + flashweave_table[otu_identifier_column].astype(str)
   my_otu_table['microbetag_id']           = flashweave_table[otu_identifier_column]
   file_to_save                            = os.path.join(FLASHWEAVE_OUTPUT_DIR, "otu_table_flashweave_format.tsv")

   flashweave_table.to_csv(file_to_save, sep ='\t', index = False)

   return my_otu_table


def edge_list_of_ncbi_ids(edgelist, species_to_ncbi_ids):

   f = open(edgelist, "r")
   associations = f.readlines()

   counter = 0

   species_present_ncbi_ids_as_dict = species_to_ncbi_ids.to_dict(orient = 'records')

   associated_pairs = {}

   for association in associations[2:]:

      taxon_a = association.split("\t")[0]
      taxon_b = association.split("\t")[1]

      taxon_a_index = None
      taxon_b_index = None

      for index, entry in enumerate(species_present_ncbi_ids_as_dict):

         if taxon_a == entry['microbetag_id']:
            taxon_a_index = index

         if taxon_b == entry['microbetag_id']:
            taxon_b_index = index

      if taxon_a_index != None and taxon_b_index != None: 
         associated_pairs[str(len(associated_pairs))] = {}
         associated_pairs[str(len(associated_pairs) - 1)]['taxon_1'] = {}
         associated_pairs[str(len(associated_pairs) - 1)]['taxon_1'] = species_present_ncbi_ids_as_dict[taxon_a_index]
         associated_pairs[str(len(associated_pairs) - 1)]['taxon_2'] = {}
         associated_pairs[str(len(associated_pairs) - 1)]['taxon_2'] = species_present_ncbi_ids_as_dict[taxon_b_index]


   return associated_pairs




# number_of_commented_lines = count_comment_lines(my_otu_table, my_taxonomy_column)

# otu_table = pd.read_csv(my_otu_table, sep = "\t", skiprows= number_of_commented_lines)



# if species_to_ncbi_ids.isin([taxon_a]).any().any() and species_to_ncbi_ids.isin([taxon_b]).any().any():
   # df2_a = species_to_ncbi_ids.loc[species_to_ncbi_ids['microbetag_id'] == taxon_a]
   # df2_b = species_to_ncbi_ids.loc[species_to_ncbi_ids['microbetag_id'] == taxon_b]
   # df2   = pd.concat([df2_a.reset_index(), df2_b.reset_index()], axis = 1)

