#!/usr/bin/bash

level=7
levelName="species"

awk -F" " '{print $1}' qiime2_consensus_taxonomy_7_levels.tsv > accessions

awk -v var=$level -F";" '{print $var}' qiime2_consensus_taxonomy_7_levels.tsv > taxa2map

paste -d "\t"  accessions taxa2map > access2taxa_silva132.tsv

../map_accession_to_ncbi_id.awk taxmap_embl_ssu_ref_132.txt access2taxa_silva132.tsv >  taxa2ncbiId.tsv.tmp

grep -v -e 'unculture\|uncultivated\|metagenome\|Ambiguous_taxa' taxa2ncbiId.tsv.tmp > taxa_names_to_ncbi_id.tsv

more taxa_names_to_ncbi_id.tsv | sort | uniq | sort > taxa_names_to_ncbi_id.tsv.uniq

mv taxa_names_to_ncbi_id.tsv.uniq qiime2_$levelName\_names_to_ncbi_id.tsv 

rm taxa2ncbiId.tsv.tmp taxa2map accessions  access2taxa_silva132.tsv taxa_names_to_ncbi_id.tsv 


