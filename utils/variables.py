"""
Global variables
"""

__author__  = 'Haris Zafeiropoulos'
__email__   = 'haris-zaf@hcmr.gr'
__status__  = 'Development'
__license__ = 'GPLv3'


import os
import sys 
import yaml
from .arguments import *


# Read configuration file
with open(args.conf, "r") as ymlfile:
   cfg = yaml.load(ymlfile, Loader = yaml.FullLoader)


BASE      = os.getcwd() 

if cfg['approach']['container']:
   IO_PATH   = "/mnt/"

else:
   IO_PATH   = cfg['approach']['io_path']


OTU_TABLE = os.path.join(IO_PATH, cfg['otu_table'])

OUT_DIR   = os.path.join(IO_PATH, cfg['output_directory'])


TAX_COL   = cfg['taxonomy_column_name']
OTU_COL   = cfg['otu_identifier_column']
COM_CHAR  = cfg['comments_character']

if cfg['column_names_are_in'] == True: 
   COM_HEAD  = '"' + 'last_comment_line' + '"'


EDGE_LIST = cfg['edge_list']

METADATA_FILE = cfg['metadata_file']




# Paths
REF_DBS                     = os.path.join(BASE, "ref-dbs")
SILVA_SPECIES_NCBI_ID       = os.path.join(REF_DBS, "silva/species_names_to_ncbi_id.tsv")
KMODULES_DEFINITIONS        = os.path.join(REF_DBS, "kegg_mappings/module_definitions.tsv")
KMODULES_DEFINITIONS_PARSED = os.path.join(REF_DBS, "kegg_mappings/module_definition_map.json")
TOOLS                       = os.path.join(BASE, "tools")

# FlashWeave
FLASHWEAVE_SCRIPT           = os.path.join(TOOLS, "flashweave/flashweave.jl")
FLASHWEAVE_OUTPUT_DIR       = os.path.join(OUT_DIR, "flashweave")
FLASHWEAVE_TMP_INPUT        = os.path.join(FLASHWEAVE_OUTPUT_DIR, "otu_table_flashweave_format.tsv")
FLASHWEAVE_EDGELIST         = os.path.join(FLASHWEAVE_OUTPUT_DIR, "network_output.edgelist")

# FAPROTAX
FAPROTAX_SCRIPT             = os.path.join(TOOLS, "faprotax/collapse_table.py")
FAPROTAX_DB                 = os.path.join(TOOLS, "faprotax/FAPROTAX.txt")

FAPROTAX_OUTPUT_DIR         = os.path.join(OUT_DIR, "faprotax")
FAPROTAX_FUNCT_TABLE        = os.path.join(FAPROTAX_OUTPUT_DIR, "functional_otu_table.tsv")
FAPROTAX_SUB_TABLES         = os.path.join(FAPROTAX_OUTPUT_DIR, "sub_tables")

# BugBase
BUGBASE_OUTPUT              = os.path.join(OUT_DIR, "bugbase")
BUGBASE_OPTIONAL            = cfg['bugbase_opt']




