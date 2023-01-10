"""
Parse user's input OTU table
"""
import pandas as pd
import logging
import os, sys
from utils.variables import * 
import networkx as nx

def count_comment_lines(my_otu_table, my_taxonomy_column):
   """
   Get the number of rows of the OTU table that should be skipped 
   """
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

def map_otu_to_ncbi_tax_level_and_id(otu_table, my_taxonomy_column, otu_identifier_column):
   """ 
   Parse user's OTU table and the Silva database to get to add 2 extra columns in the OTU table: 
   1. the lowest taxonomic level of the taxonomy assigned in an OTU, for which an NCBI Taxonomy id exists (e.g., "genus")
   2. the corresponding NCBI Taxonomy Id (e.g., "343")
   """
   taxonomies = otu_table.filter([otu_identifier_column, my_taxonomy_column])

   # Split the taxonomy column and split it based on semicolumn! 
   splitted_taxonomies         = taxonomies[my_taxonomy_column].str.split(TAX_DELIM, expand = True)
   splitted_taxonomies.columns = ["Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
   splitted_taxonomies[otu_identifier_column] = taxonomies[otu_identifier_column]

   # Make sure there is no white space character in the Species, Genus and Family levels
   splitted_taxonomies['Species'] = splitted_taxonomies['Species'].str.strip()
   splitted_taxonomies['Genus']   = splitted_taxonomies['Genus'].str.strip()
   splitted_taxonomies['Family']  = splitted_taxonomies['Family'].str.strip()

   # Read the species, genus and family files of Silva db with their NCBI IDs
   silva_species_ncbi_id                   = pd.read_csv(SILVA_SPECIES_NCBI_ID, sep = "\t")
   silva_species_ncbi_id.columns = ['Species', 'ncbi_tax_id']
   silva_species_ncbi_id['Species'].str.strip()

   silva_genus_ncbi_id                   = pd.read_csv(SILVA_GENUS_NCBI_ID, sep = "\t")
   silva_genus_ncbi_id.columns = ['Genus', 'ncbi_tax_id']
   silva_genus_ncbi_id['Genus'].str.strip()

   silva_family_ncbi_id                   = pd.read_csv(SILVA_FAMILY_NCBI_ID, sep = "\t")
   silva_family_ncbi_id.columns = ['Family', 'ncbi_tax_id']
   silva_family_ncbi_id['Family'].str.strip()

   # Build a dataframe for the Species, Genus and Family taxonomies present on the OTU table along with their corresponding NCBI Tax IDs
   species_present  = silva_species_ncbi_id.merge(splitted_taxonomies, on = ['Species'])
   species_present  = species_present[['Species', 'ncbi_tax_id', otu_identifier_column]]
   species_present.insert(0,"ncbi_tax_level", "species") 

   genera_present    = silva_genus_ncbi_id.merge(splitted_taxonomies, on = ['Genus'])
   genera_present    = genera_present[['Genus', 'ncbi_tax_id', otu_identifier_column]]
   genera_present.insert(0,"ncbi_tax_level", "genus") 
   
   families_present = silva_family_ncbi_id.merge(splitted_taxonomies, on = ['Family'])
   families_present = families_present[['Family', 'ncbi_tax_id', otu_identifier_column]]
   families_present.insert(0,"ncbi_tax_level", "family") 

   # Buid a merged df with the NCBI level and the corresponding tax id
   merged_taxonomies = splitted_taxonomies.set_index(OTU_COL)
   merged_taxonomies.insert(0,"ncbi_tax_id", "NaN")
   merged_taxonomies.insert(0,"ncbi_tax_level", "NaN")

   species_present = splitted_taxonomies.set_index(OTU_COL)
   genera_present = genera_present.set_index(OTU_COL)
   families_present = families_present.set_index(OTU_COL)

   merged_taxonomies.update(families_present)
   merged_taxonomies.reset_index(inplace=True)

   merged_taxonomies.update(genera_present)
   merged_taxonomies.reset_index(inplace=True)

   merged_taxonomies.update(species_present)
   merged_taxonomies.reset_index(inplace=True)
   merged_taxonomies.drop(columns=["level_0", "index"],  inplace=True)
 
   merged_taxonomies[OTU_COL] = merged_taxonomies[OTU_COL].apply(str)
 
   return merged_taxonomies

def otu_faprotax_functions_assignment(path_to_subtables):
   """
   Parse the sub tables of the faprotax analysis 
   to assign the biological processes related to each OTU 
   """
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

def build_annotated_graph(edgelist, otu2ncbi, **kwargs):

   G = nx.Graph()

   nodes = []
   edges = []

   for otu_id, ncbi_elements in otu2ncbi.items(): 

      otu_id = str(otu_id)

      node               = {}
      node["data"]       = {}
      node["data"]["id"] = otu_id
      try:
         node["data"]["ncbi-id"] = ncbi_elements[otu_id]["ncbi_tax_id"]
      except:
         node["data"]["ncbi-id"] = "NA"

      for resource, annotations in kwargs.items():
         if resource == "fapr": 
            try: 
               node["data"]["faprotax-annotations"] = annotations[otu_id]
            except:
               node["data"]["faprotax-annotations"] = "None"
         """ 
         when more resources are available, just expand this try-except case. 
         """
      nodes.append(node)

   for association, edge_elements in edgelist.items():

      edge               = {} 
      edge["data"]       = {}
      edge["data"]["id"] = association
      edge["data"]["source"] = edge_elements["taxon_A"]
      edge["data"]["target"] = edge_elements["taxon_B"]
      if not EDGE_LIST: 
         edge["data"]["FlashWeave-weight"] = edge_elements["evidence"]

def is_tab_separated(my_otu_table, my_taxonomy_column):
   """
   Read the OTU table and make sure it is tab separated and not empty
   """

   number_of_commented_lines = count_comment_lines(my_otu_table, my_taxonomy_column)

   try:
      otu_table = pd.read_csv(my_otu_table, sep = OTU_TABLE_DELIM, skiprows= number_of_commented_lines)
   except:
      logging.error("The OTU table provided is not a tab separated file. Please convert your OTU table to .tsv or .csv format.")

   if otu_table.shape[1] < 2 :
      logging.error("The OTU table you provide has no records.")

   return otu_table

def ensure_flashweave_format(my_otu_table, my_taxonomy_column, otu_identifier_column):
   """
   Build an OTU table that will be in a FlashWeave-based format. 
   """
   flashweave_table = my_otu_table.drop(my_taxonomy_column, axis = 1)
   float_col                   = flashweave_table.select_dtypes(include=['float64']) 
   
   for col in float_col.columns.values:
      flashweave_table[col] = flashweave_table[col].astype('int64')

   flashweave_table[otu_identifier_column] = 'microbetag_' + flashweave_table[otu_identifier_column].astype(str)
   my_otu_table['microbetag_id']                     = flashweave_table[otu_identifier_column]

   file_to_save  = os.path.join(FLASHWEAVE_OUTPUT_DIR, "otu_table_flashweave_format.tsv")
   flashweave_table.to_csv(file_to_save, sep ='\t', index = False)

   return my_otu_table

def edge_list_of_ncbi_ids(edgelist, otu_table_with_ncbi_ids):
   """
   Read an edge list and build a dataframe with the corresponding ncbi ids for each pair 
   if and only if, both OTUs have been mapped to a NCBI tax id 
   e.g.
       ncbi_tax_id_a       ncbi_tax_level_a     ncbi_tax_id_b      ncbi_tax_level_b
         0             838          genus                      171552           family
         1        186803     family                      186807           family
   """
   f = open(edgelist, "r")
   associations = f.readlines()
   associated_pairs = pd.DataFrame(columns=["ncbi_tax_id_a", 
                                                                                             "ncbi_tax_level_a", 
                                                                                             "ncbi_tax_id_b", 
                                                                                             "ncbi_tax_level_b"]
                                                                        )
   counter = 0
   for association in associations[2:]:

      otu_a = association.split("\t")[0]
      otu_b = association.split("\t")[1]

      if not EDGE_LIST:
         otu_a = otu_a[11:]
         otu_b = otu_b[11:]

      ncbi_id_otu_a      = otu_table_with_ncbi_ids.loc[otu_table_with_ncbi_ids[OTU_COL] == otu_a]["ncbi_tax_id"]
      ncbi_level_otu_a = otu_table_with_ncbi_ids.loc[otu_table_with_ncbi_ids[OTU_COL] == otu_a]["ncbi_tax_level"]
      ncbi_id_otu_b      = otu_table_with_ncbi_ids.loc[otu_table_with_ncbi_ids[OTU_COL] == otu_b]["ncbi_tax_id"]
      ncbi_level_otu_b = otu_table_with_ncbi_ids.loc[otu_table_with_ncbi_ids[OTU_COL] == otu_b]["ncbi_tax_level"]

      if ncbi_id_otu_a.item() != "NaN" and ncbi_id_otu_b.item() !="NaN":
         df = pd.DataFrame(
            { "ncbi_tax_id_a": ncbi_id_otu_a.item(),
              "ncbi_tax_level_a": ncbi_level_otu_a.item(),
              "ncbi_tax_id_b": ncbi_id_otu_b.item(),
              "ncbi_tax_level_b": ncbi_level_otu_b.item()
            }, 
            index=[counter]
         )
         associated_pairs = pd.concat([associated_pairs, df], axis=0)
         counter += 1

   # In the otu table pd we have floats and strings in the Id column, so we need to make sure we now have only strings 
   associated_pairs["ncbi_tax_id_a"] = associated_pairs["ncbi_tax_id_a"].apply(int)
   associated_pairs["ncbi_tax_id_a"] = associated_pairs["ncbi_tax_id_a"].apply(str)
   associated_pairs["ncbi_tax_id_b"] = associated_pairs["ncbi_tax_id_b"].apply(int)
   associated_pairs["ncbi_tax_id_b"] = associated_pairs["ncbi_tax_id_b"].apply(str)

   return associated_pairs

def get_species(my_otu_table, my_taxonomy_column, otu_identifier_column):
   """
   DEPRECATED FUNCTION
   """
   try: 

      """
      pd.filter(): 
         Subset rows or columns of dataframe according to labels in the specified index. 
         Note that this routine does not filter a dataframe on its contents. 
         The filter is applied to the labels of the index.
      """
      # Keep only otu ids and taxonomies assigned
      otu_id_and_taxonomy = my_otu_table.filter([otu_identifier_column, my_taxonomy_column, 'microbetag_id']) 

      # Keep only the last level of lineage assigned and remove white spaces if any before the species name assigned
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

def edgelist_to_ncbi_ids(edgelist, otu_2_tax_level):
   """
   DEPRECATED FUNCTION
   """

   """
   Flag associations of edge list based on whether both taxa are at the species/strain level 
   or at a higher one.
   """
   if EDGE_LIST: 
      associations = pd.read_csv(edgelist, sep = "\t", skiprows = 2)
   else: 
      associations = pd.read_csv(edgelist, sep = "\t", skiprows = 2)

   associations.columns = ['taxon_A', 'taxon_B', 'evidence']

   associations_dict = associations.to_dict('index')

   for k,v in associations_dict.items(): 

      otu_1 = v['taxon_A']
      otu_2 = v['taxon_B']

      try:
         if otu_2_tax_level[otu_1]['tax_level'] == "Species" and otu_2_tax_level[otu_2]['tax_level'] == "Species":
            associations_dict[k]['tax_level'] = "species"
         else:
            associations_dict[k]['tax_level'] = "genus or family"
      except:
         associations_dict[k]['tax_level'] = "NA"

   return associations_dict
