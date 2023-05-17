"""
Global variables
"""
import os
import sys 
import yaml
from .arguments import *


BASE = "/".join(os.path.abspath(__file__).split("/")[:-2])

# DATABASE 
HOST      = "localhost"
USER_NAME = os.environ.get('DB_USER')
PASSWORD  = os.environ.get('DB_PASSWORD')
DB_NAME   = 'microbetagDB'

# Paths to reference database and mapping folders 
REF_DBS  = os.path.join(BASE, "ref-dbs")
MAPPINGS = os.path.join(BASE, "mappings")

# GTDB 
GTDB_METADATA = os.path.join(MAPPINGS, "GTDB_QUALITY_REPRESENTATIVE_GENOMES_v202")

# NCBI Ids paths
SPECIES_NCBI_IDS  = os.path.join(MAPPINGS, "species2ncbiId.tsv")
GENERA_NCBI_IDS   = os.path.join(MAPPINGS, "genera2ncbiId.tsv")
FAMILIES_NCBI_IDS = os.path.join(MAPPINGS, "families2ncbiId.tsv")

# KEGG paths - Pathway complmementarity module
KMODULES_DEFINITIONS        = os.path.join(REF_DBS, "kegg_mappings/module_definitions.tsv")
KMODULES_DEFINITIONS_PARSED = os.path.join(REF_DBS, "kegg_mappings/module_definition_map.json")
ALL_GENOMES_MODULES         = os.path.join(REF_DBS, "all_genomes_modules")

# Tools
TOOLS             = os.path.join(BASE, "tools")
FLASHWEAVE_SCRIPT = os.path.join(TOOLS, "flashweave/flashweave.jl")
FAPROTAX_SCRIPT   = os.path.join(TOOLS, "faprotax/collapse_table.py")
FAPROTAX_DB       = os.path.join(TOOLS, "faprotax/FAPROTAX.txt")



# Read configuration file
check = True
try:
   with open(args.conf, "r") as ymlfile:
      cfg = yaml.load(ymlfile, Loader = yaml.FullLoader)
except:
   check = False

if check:

   if cfg["approach"]["container"]:
      IO_PATH   = "/mnt/"
   else:
      IO_PATH = cfg["approach"]["io_path"]


   if cfg[ "taxonomy" ] == "GTDB":
      GTDB_ACCESSION_NCBI_TAX_IDS = os.path.join(MAPPINGS, "species2ncbiId2accession.tsv")
   elif cfg["taxonomy"] == "dada2":
      GTDB_ACCESSION_NCBI_TAX_IDS = os.path.join(MAPPINGS, "dada2ncbi2accession.tsv")
   elif cfg["taxonomy"] == "qiime2":
      GTDB_ACCESSION_NCBI_TAX_IDS = os.path.join(MAPPINGS, "qiime2species2ncbi2accession.tsv")
   else:
      GTDB_ACCESSION_NCBI_TAX_IDS = os.path.join(MAPPINGS, "gtdbSpecies2ncbiId2accession.tsv")



   OTU_TABLE       = os.path.join(IO_PATH, cfg["otu_table"])
   OTU_TABLE_DELIM = cfg["otu_table_delim"]
   OUT_DIR         = os.path.join(IO_PATH, cfg["output_directory"])

   TAX_COL   = cfg["taxonomy_column_name"]
   TAX_DELIM = cfg["taxonomy_delimeter"]
   OTU_COL   = cfg["otu_identifier_column"]
   COM_CHAR  = cfg["comments_character"]

   if cfg["column_names_are_in"] == True: 
      COM_HEAD  = '"' + "last_comment_line" + '"'
   else:
      COM_HEAD = ""

   EDGE_LIST     = os.path.join(IO_PATH, cfg[ "edge_list" ]) if cfg[ "edge_list" ] != None else False
   METADATA_FILE = cfg[ "metadata_file" ]

   # Steps 
   PATHWAY_COMPLEMENTARITY = True if cfg[ "pathway_complementarity" ] != False else False
   PHEN_DB                 = True if cfg[ "PhenDB" ]   != False else False
   BUGBASE                 = True if cfg[ "BugBase" ]  != False else False
   FAPROTAX                = True if cfg[ "FAPROTAX" ] != False else False


   USERS_TAXONOMY          = cfg[ "taxonomy" ]


   # FlashWeave
   FLASHWEAVE_OUTPUT_DIR       = os.path.join(OUT_DIR, "flashweave")
   FLASHWEAVE_TMP_INPUT        = os.path.join(FLASHWEAVE_OUTPUT_DIR, "otu_table_flashweave_format.tsv")
   FLASHWEAVE_EDGELIST         = os.path.join(FLASHWEAVE_OUTPUT_DIR, "network_output.edgelist")

   # FAPROTAX
   FAPROTAX_OUTPUT_DIR         = os.path.join(OUT_DIR, "faprotax")
   FAPROTAX_FUNCT_TABLE        = os.path.join(FAPROTAX_OUTPUT_DIR, "functional_otu_table.tsv")
   FAPROTAX_SUB_TABLES         = os.path.join(FAPROTAX_OUTPUT_DIR, "sub_tables")

   # BugBase
   BUGBASE_PATH     = os.getenv("BUGBASE_PATH")
   BUGBASE_OPTIONAL = cfg[ "bugbase_opt" ]
   BUGBASE_SCRIPT   = os.path.join(BUGBASE_PATH, "bin/run.bugbase.r")
   BUGBASE_TMP      = os.path.join(OUT_DIR, "tmp_bugbase_otu_table.txt")
   BUGBASE_OUTPUT   = os.path.join(OUT_DIR, "bugbase")
 
   # PhenDB
   PHEN_OUTPUT_DIR = os.path.join(OUT_DIR, "phen_predictions")

