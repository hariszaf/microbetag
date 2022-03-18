"""
Global variables
"""
import os
import sys # even if not needed here, it is used on the microbetag.py script

# A relative import specifies the resource to be imported relative to the current location; 
# the location where the import statement is.
# A single dot means that the module or package referenced is in the same directory as the current location. 
# Two dots mean that it is in the parent directory of the current location—that is, the directory above. 
# Three dots mean that it is in the grandparent directory, and so on.
# For more: 
# # For more have a look at: https://realpython.com/absolute-vs-relative-python-imports/#syntax-and-practical-examples_1
# Remove . if you want to test the script on its own - not calling it from the microbetag.py script
from .arguments import *

BASE      = os.getcwd() 
OUT_DIR   = args.o

EDGE_LIST = args.el

OTU_TABLE = args.i
TAX_COL   = args.t
OTU_COL   = args.c
COM_CHAR  = args.com
if args.ch != False: 
   COM_HEAD  = '"' + 'last_comment_line' + '"'


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
         if TAX_COL is None:
            TAX_COL = line[-1].rstrip()
            # TAX_COL = '"' + line[-1].rstrip() + '"'
         if OTU_COL is None: 
            OTU_COL = line[0].rstrip()
         break
      if args.com is None:
         COM_CHAR = '"' + line[0][0] + '"'


REF_DBS                     = os.path.join(BASE, "ref-dbs")
SILVA_SPECIES_NCBI_ID       = os.path.join(REF_DBS, "silva/species_names_to_ncbi_id.tsv")
KMODULES_DEFINITIONS        = os.path.join(REF_DBS, "kegg_mappings/module_definitions.tsv")
KMODULES_DEFINITIONS_PARSED = os.path.join(REF_DBS, "kegg_mappings/module_definition_map.json")


TOOLS                       = os.path.join(BASE, "tools")

FLASHWEAVE_SCRIPT           = os.path.join(TOOLS, "flashweave/flashweave.jl")
FLASHWEAVE_OUTPUT_DIR       = os.path.join(OUT_DIR, "flashweave")
FLASHWEAVE_TMP_INPUT        = os.path.join(FLASHWEAVE_OUTPUT_DIR, "otu_table_flashweave_format.tsv")
# print("#@#$@#$ -->", FLASHWEAVE_TMP_INPUT)

FAPROTAX_SCRIPT             = os.path.join(TOOLS, "faprotax/collapse_table.py")
FAPROTAX_DB                 = os.path.join(TOOLS, "faprotax/FAPROTAX.txt")







FAPROTAX_OUTPUT_DIR         = os.path.join(OUT_DIR, "faprotax")
FAPROTAX_FUNCT_TABLE        = os.path.join(FAPROTAX_OUTPUT_DIR, "functional_otu_table.tsv")
FAPROTAX_SUB_TABLES         = os.path.join(FAPROTAX_OUTPUT_DIR, "sub_tables")



