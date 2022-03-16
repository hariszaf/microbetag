#!/usr/bin/env python

"""
Global variables
"""
import os


BASE                        = os.path.dirname(os.getcwd())
REF_DBS                     = os.path.join(BASE, "ref-dbs")
SILVA_SPECIES_NCBI_ID       = os.path.join(REF_DBS, "silva/species_names_to_ncbi_id.tsv")
KMODULES_DEFINITIONS        = os.path.join(REF_DBS, "kegg_mappings/module_definitions.tsv")
KMODULES_DEFINITIONS_PARSED = os.path.join(REF_DBS, "kegg_mappings/module_definition_map.json")




