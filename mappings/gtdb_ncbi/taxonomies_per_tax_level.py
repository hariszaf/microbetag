#!/usr/bin/env python3

__name__    = "microbetag"
__author__  = "Haris Zafeiropoulos"
__email__   = "haris.zafeiropoulos@kuleuven.be"
__license__ = "GPLv3"
__version__ = "v0.1.0"

"""
Aim: 
This script gets as input the metadata files describing GTDB (e.g. ar122_metadata_r202.tar.gz)
one can find from the GTDB Downloads (https://data.gtdb.ecogenomic.org/releases/release202/202.0/)
and returns the unique taxa names included on the database per taxonomic level (i.e., phylum, class ets.)

Under its current version, microbetag makes use of GTDB v202.

Usage: ./taxonomies_per_tax_level.py <metadata_file> <output_prefix> <column_number_with_GTDB_taxonomy>

"""

import sys, os

metadata_file = open(sys.argv[1], "r")
output_prefix = sys.argv[2]
col_number    = int(sys.argv[3]) - 1

phyla    = set()
classes  = set()
orders   = set()
families = set()
genera   = set()
species  = []

species_accession = {}


counter = 0
for line in metadata_file: 

    if counter == 0:
        counter += 1
        continue

    line            = line.split("\t")
    gtdb_taxonomy   = line[ col_number ].split(";")

    phyla.add(gtdb_taxonomy[1])
    classes.add(gtdb_taxonomy[2])
    orders.add(gtdb_taxonomy[3])
    families.add(gtdb_taxonomy[4])
    genera.add(gtdb_taxonomy[5])
    
    if gtdb_taxonomy[6] not in species:
        species.append(gtdb_taxonomy[6])
        species_accession[ gtdb_taxonomy[6] ] = line[0]

#with open( "".join( [output_prefix, "_phyla.tsv"] ), "w" ) as f:
#    for entry in phyla:
#        f.write(entry + "\n")
#
#with open( "".join( [output_prefix, "_classes.tsv"] ), "w" ) as f:
#    for entry in classes:
#        f.write(entry + "\n")
#
#with open( "".join( [output_prefix, "_orders.tsv"] ), "w") as f:
#    for entry in orders:
#        f.write(entry + "\n")
#
#with open( "".join( [output_prefix, "_families.tsv"] ), "w") as f:
#    for entry in families:
#        f.write(entry + "\n")
#
#with open( "".join( [output_prefix, "_genera.tsv"] ), "w") as f:
#    for entry in genera:
#        f.write(entry + "\n")
#
with open( "".join( [output_prefix, "_species.tsv"]), "w") as f:
    for entry in species:
        f.write(species_accession[entry] + "\t" + entry + "\n")
