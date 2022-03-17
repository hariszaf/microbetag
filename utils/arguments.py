"""
Parse arguments
"""

import argparse

parser = argparse.ArgumentParser()

# parser.add_argument('-conf', '--configuration_file', metavar='configuration_file', dest='conf',
#                     required = True,
#                     help='YAML configuration file with the parameter values')

# [TO_REMEMBER!] Replace required as True when ready!
parser.add_argument('-i', '--otu_table', metavar = 'input_otu_table', dest = 'i',
                    type = str, required = True,
                    help = 'OTU table with taxonomies assigned using the Silva database. In case other database was used, the OTU table needs to have an extra column specifying the NCBI Taxonomy id of the assignments that are up to the species level')

# If none, the last column of the OTU table will be used as such 
parser.add_argument('-t', '--taxonomy_column_name', metavar='taxonomy_column', dest='t',
                    type=str, required = False,
                    help='exact name of the column denoting the taxonomy assignment in the OTU table')

# If none, the first column of the OTU table will be used as such 
parser.add_argument('-c', '--otu_identifier_column', metavar='taxonomy_column', dest='c',
                    type=str, required = False,
                    help='exact name of the column denoting the taxonomy assignment in the OTU table')

# If none, the first character of the last line of the non taxa lines will be used
parser.add_argument('-com', '--comments_character', metavar='comments_character', dest='com',
                    type=str, required = False,
                    help='the character denoting commented lines on the OTU table')

# Metadata table as described in the BugBase documentation
parser.add_argument('-m', '--mapping_file', metavar = 'mapping_file', dest = 'm',
                    type = str, required = False,
                    help = 'Metadata table as described in the BugBase documentation. See: https://bugbase.cs.umn.edu/documentation.html')


# If none, the output directory will be built in the current working directory
parser.add_argument('-o', '--outdir', metavar='output_dir', dest='o',
                    type = str, required = False, default = "output",
                    help = 'output directory')

# By default, microbetag assumes the OTU table is in Silva v.132
parser.add_argument('-is', '--is_silva', metavar='is_silva_taxonomy', dest='is',
                    type = str, required = False, default = "True",
                    help = 'threads that will be used')

args = parser.parse_args()
