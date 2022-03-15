#!/usr/bin/env python 


"""
Function to read an OTU table and assign 
an NCBI Taxonomy Id to each assigned taxonomy
"""


esearch -db taxonomy -query "enterococcus" | esummary | xtract -pattern DocumentSummary -element ScientificName,Rank,TaxId



"""
Function to add a NCBI Taxonomy Id to the species 
of an OTU table.
"""
def read_otu_table(otu_table, taxonomy_scheme):


   # BIOM VS TSV

   # SILVA multiple versions 


   for 