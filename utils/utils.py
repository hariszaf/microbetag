"""
Parse user's input OTU table
"""
import pandas as pd
import numpy as np
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

   [TODO] The code of this function can be more clear by writing mspecies, species and genera levels as the on of the family.
   """
   taxonomies = otu_table.filter([otu_identifier_column, 
                                  my_taxonomy_column, 
                                  "microbetag_id"])

   # Split the taxonomy column and split it based on semicolumn! 
   splitted_taxonomies         = taxonomies[my_taxonomy_column].str.split(TAX_DELIM, expand = True)
   splitted_taxonomies.columns = ["Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
   splitted_taxonomies[otu_identifier_column] = taxonomies[otu_identifier_column]
   splitted_taxonomies["microbetag_id"] = taxonomies["microbetag_id"]

   splitted_taxonomies_c = splitted_taxonomies.copy()


   # Read the species, genus and family files of Silva db with their NCBI IDs
   """
   The SPECIES_NCBI_ID, GENERA_NCBI_IDS and FAMILIES_NCBI_IDS files are like this: 
   Species  ncbi_tax_id
   D_6__Abiotrophia sp. oral clone OH2A       319434
   """
   gtdb_accession_ids         = pd.read_csv(GTDB_ACCESSION_NCBI_TAX_IDS, sep = "\t")
   gtdb_accession_ids.columns = ["Species", "ncbi_tax_id", "gtdb_gen_repr"]
   gtdb_accession_ids["Species"].str.strip()

   species_ncbi_id         = pd.read_csv(SPECIES_NCBI_IDS, sep = "\t")
   species_ncbi_id.columns = ['Species', 'ncbi_tax_id']
   species_ncbi_id['Species'].str.strip()

   genera_ncbi_id         = pd.read_csv(GENERA_NCBI_IDS, sep = "\t")
   genera_ncbi_id.columns = ['Genus', 'ncbi_tax_id']
   genera_ncbi_id['Genus'].str.strip()

   family_ncbi_id         = pd.read_csv(FAMILIES_NCBI_IDS, sep = "\t")
   family_ncbi_id.columns = ['Family', 'ncbi_tax_id']
   family_ncbi_id['Family'].str.strip()

   # Build a dataframe for the Species, Genus and Family taxonomies present 
   # on the OTU table along with their corresponding NCBI Tax IDs

   # GTDB accession ids
   genomes_present = gtdb_accession_ids.merge(
                                             splitted_taxonomies, 
                                             on = ["Species"]
                                             )[["ncbi_tax_id","gtdb_gen_repr", otu_identifier_column]]

   splitted_taxonomies = pd.merge(genomes_present, splitted_taxonomies, 
                                  on = otu_identifier_column, how = 'outer')
   splitted_taxonomies.loc[ splitted_taxonomies["ncbi_tax_id"].notnull(), "ncbi_tax_level" ] = "mspecies"

   mspecies = splitted_taxonomies[[otu_identifier_column, "gtdb_gen_repr"]]

   # Species
   species_present  = species_ncbi_id.merge(splitted_taxonomies.query('ncbi_tax_level != "mspecies"'), 
                                                                     on = ['Species'], 
                                                                     suffixes = ('_species', '_mspecies')
                                          )[
                                             [otu_identifier_column, "ncbi_tax_id_species","ncbi_tax_level"]
                                          ]

   species_present.loc[ species_present["ncbi_tax_level"] != "mspecies", "ncbi_tax_level"] = "species"

   splitted_taxonomies = pd.merge(splitted_taxonomies, species_present,
                                 on  = otu_identifier_column, 
                                 how ="outer", 
                                 suffixes = ('_species', '_mspecies')
                                 )

   splitted_taxonomies["ncbi_tax_id"]    = splitted_taxonomies["ncbi_tax_id"].fillna(splitted_taxonomies["ncbi_tax_id_species"])
   splitted_taxonomies["ncbi_tax_level"] = splitted_taxonomies["ncbi_tax_level_mspecies"].fillna(splitted_taxonomies["ncbi_tax_level_species"])

   splitted_taxonomies = splitted_taxonomies.drop(["ncbi_tax_id_species", "ncbi_tax_level_species" , "ncbi_tax_level_mspecies"], axis = 1)

   """ 
   REMEMBER! 
   There is no one-to-one relationship between a NCBI Taxonomy Id and a representative genome 

   We do not have a hit for ALL NCBI TAX IDs from the species found in our experiment, aka there s 
   no GTDB representative genome for all NCBI Taxonomy Ids, e.g.  Bradyrhizobium sp. J81, NCBI Tax Id: 656743 
   """

   # Genera
   pd_genera         = splitted_taxonomies[[otu_identifier_column, "ncbi_tax_level", "ncbi_tax_id", "gtdb_gen_repr", "Genus"]]
   genera_present    = genera_ncbi_id.merge(pd_genera, on = ['Genus'], suffixes=("_gen", "_over"), how="right")

   genera_present.loc[ genera_present["ncbi_tax_id_gen"].notnull(), "ncbi_tax_level_gen" ] = "genus"

   genera_present["ncbi_tax_id"]    = genera_present['ncbi_tax_id_over'].combine_first(genera_present['ncbi_tax_id_gen'])
   genera_present["ncbi_tax_level"] = genera_present['ncbi_tax_level'].combine_first(genera_present['ncbi_tax_level_gen'])

   genera_present = genera_present.drop(["ncbi_tax_id_gen", "ncbi_tax_id_over", "Genus"], axis=1)


   # Families
   pd_families = splitted_taxonomies[[otu_identifier_column, "Family"]]
   families_present = family_ncbi_id.merge(pd_families, on = ['Family'], how="right")
   families_present.loc[ families_present["ncbi_tax_id"].notnull(), "ncbi_tax_level" ] = "family"

   families_present["ncbi_tax_id"]    = genera_present["ncbi_tax_id"].combine_first(families_present["ncbi_tax_id"])
   families_present["ncbi_tax_level"] = genera_present["ncbi_tax_level"].combine_first(families_present["ncbi_tax_level"])
   families_present = families_present.drop(["Family"], axis=1)

   # Build a unified data frame for all levels and the accession ids when available
   otu_taxid_level_repr_genome = pd.merge(splitted_taxonomies_c, families_present, on = otu_identifier_column)
   otu_taxid_level_repr_genome = pd.merge(otu_taxid_level_repr_genome, mspecies, on = otu_identifier_column, how="outer")

   repr_genomes_present = list(mspecies["gtdb_gen_repr"].dropna())

   return otu_taxid_level_repr_genome, repr_genomes_present

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

   """
   [TODO] If we come up with a microbetag database, here we need to assign to the "microbetag_id" variable
   the microbetag db ids to the corresponding taxa BUT apparently, this we ll have to be moved AFTER getting the NCBI ids.
   """
   flashweave_table[otu_identifier_column] = 'microbetag_' + flashweave_table[otu_identifier_column].astype(str)
   my_otu_table['microbetag_id']           = flashweave_table[otu_identifier_column]

   file_to_save  = os.path.join(FLASHWEAVE_OUTPUT_DIR, "otu_table_flashweave_format.tsv")
   flashweave_table.to_csv(file_to_save, sep ='\t', index = False)

   return my_otu_table

def edge_list_of_ncbi_ids(edgelist, otu_table_with_ncbi_ids):
   """
   Read an edge list and build a dataframe with the corresponding ncbi ids for each pair 
   if and only if, both OTUs have been mapped to a NCBI tax id 
   e.g.
               ncbi_tax_id_a  ncbi_tax_level_a  ncbi_tax_id_b   ncbi_tax_level_b
         0        838              genus           171552          family
         1       186803           family           186807          family
   """
   pd_edgelist = pd.read_csv(edgelist, sep="\t", skiprows=2, header=None)
   pd_edgelist.columns = ["node_a", "node_b", "score"]

   pd_edgelist["joint"] = pd_edgelist['node_a'].astype(str) +":"+ pd_edgelist["node_b"]

   associated_pairs_node_a = pd.merge(pd_edgelist[["node_a", "joint"]], otu_table_with_ncbi_ids[["ncbi_tax_id", "gtdb_gen_repr", "ncbi_tax_level", "microbetag_id"]],
                                                left_on='node_a', right_on='microbetag_id', how="inner").drop(["microbetag_id"], axis=1)
   associated_pairs_node_a.rename(columns = {
      "ncbi_tax_level": "ncbi_tax_level_node_a",
      "gtdb_gen_repr": "gtdb_gen_repr_node_a",
      "ncbi_tax_id": "ncbi_tax_id_node_a"
   }, inplace = True)

   associated_pairs_node_b = pd.merge(pd_edgelist[["node_b", "joint", "score"]], otu_table_with_ncbi_ids[["ncbi_tax_id", "gtdb_gen_repr", "ncbi_tax_level", "microbetag_id"]],
                                                left_on='node_b', right_on='microbetag_id', how="inner").drop(["microbetag_id"], axis=1)
   associated_pairs_node_b.rename(columns = {
      "ncbi_tax_level": "ncbi_tax_level_node_b",
      "gtdb_gen_repr": "gtdb_gen_repr_node_b",
      "ncbi_tax_id": "ncbi_tax_id_node_b"
   }, inplace = True)

   associated_pairs = pd.merge(associated_pairs_node_a, associated_pairs_node_b, on="joint").drop(["joint"], axis=1)

   return associated_pairs

# NOT TO USE
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

   except KeyboardInterrupt:
      logging.error("The taxonomy column could not be retrieved from you OTU table. \n \
                     Check for any strange characters on your OTU table or its format.")

   # Load the Silva species names along with ther NCBI Taxonomy ids
   species_ncbi_id         = pd.read_csv(SPECIES_NCBI_IDS, sep = "\t")
   species_ncbi_id.columns = ['species_name', 'ncbi_tax_id']
   species_ncbi_id['species_name'].str.strip()
  
   # Make a column showing whether or not the species of each row is included on the Silva ref file
   otu_id_and_taxonomy['present'] = otu_id_and_taxonomy[my_taxonomy_column].isin(species_ncbi_id['species_name'])

   # We make a map (dictionary) keeping as a key the column that links the 2 datadrames; in this case the species name
   map_dict = dict(zip(species_ncbi_id['species_name'], species_ncbi_id['ncbi_tax_id']))
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
