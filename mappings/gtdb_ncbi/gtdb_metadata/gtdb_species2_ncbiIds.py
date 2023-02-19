#!/usr/bin/env python 

import pandas as pd 
import sys

species = open("gtdb_species.tsv","r")


gtdb = pd.read_csv("GTDB_QUALITY_REPRESENTATIVE_GENOMES_v202", sep="\t")

gtdb_ncbiId = gtdb["ncbi_taxid"]

ncbi_taxonomy = gtdb["ncbi_taxonomy"]
ncbi = ncbi_taxonomy.apply(lambda x: pd.Series(str(x).split(";")))
ncbi_species = ncbi[6]

gtdb_taxonomy = gtdb["gtdb_taxonomy"]
gtdb = gtdb_taxonomy.apply(lambda x: pd.Series(str(x).split(";")))
gtdb_species = gtdb[6]



df = pd.concat( [ gtdb_ncbiId, ncbi_species ], axis=1)
df.rename(columns={6:'ncbi_taxonomy'}, inplace=True)

df2 = pd.concat( [df, gtdb_species], axis = 1)
df2.rename(columns={6:'gtdb_taxonomy'}, inplace=True)


output = open("gtdb_species2ncbi_ids.tsv", "w")

for line in species:

    species_name = line[:-1]

    ncbi_match = df.loc[ df2["ncbi_taxonomy"] == species_name ]
    gtdb_match = df.loc[ df2["gtdb_taxonomy"] == species_name ]

    if not gtdb_match.empty:
        output.write(species_name + "\t" + str(gtdb_match["ncbi_taxid"].to_list()[0]) + "\n")
        continue

    if not ncbi_match.empty:
        output.write(species_name + "\t" + str(ncbi_match["ncbi_taxid"].to_list()[0]) + "\n")

