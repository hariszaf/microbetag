Based on the [Qiime files of Silva v.132](https://www.arb-silva.de/download/archive/qiime/) we built 
a two-column file (`species_names_to_ncbi_id.tsv`) with the species names of the 7-level 
taxonomy scheme and their corresponding NCBI Taxonomy id. 

 
We consider as *species names* only those entries where at the 7th level 
of the full lineage, an actual species name is available. 
For example, the species name `;D_6__Vannella sp. LITHOV` is included as 
it comes from a full lineage:  
```
D_0__Eukaryota;D_1__Amoebozoa;D_2__Discosea;D_3__Flabellinia;D_4__Van
nellida;D_5__Platyamoeba;D_6__Vannella sp. LITHOV
```

On the contrary, in the case of: 

```
AB742067.1.1456 D_0__Bacteria;D_1__Firmicutes;D_2__Negativicutes;D_3__Selenomonadales
;D_4__Veillonellaceae;D_5__uncultured;D_6__uncultured Firmicutes bacterium
```

where at the 7th level we find the term `D_6__uncultured Firmicutes bacterium`,
this is not included. 


> It is the entries at the file `species_names_to_ncbi_id.tsv` that will be considered for the pathway complementarity module

## Steps

Using Silva v132 and based on the
[Qiime files of Silva v.132](https://www.arb-silva.de/download/archive/qiime/) we ran the following commands to build our 2-columns file:


```bash=
grep -v "Eukaryota" cons_taxonomy_7_levels_silva_132.txt > cons_taxonomy_7_levels_silva_132_no_eukaryotes.txt

awk -F" " '{print $1}' cons_taxonomy_7_levels_silva_132_no_eukaryotes.txt > accessions

awk -F";" '{print $NF}' cons_taxonomy_7_levels_silva_132_no_eukaryotes.txt > species

paste -d "\t"  accessions species > access2species_silva132.tsv

./map_accession_to_ncbi_id.awk taxmap_embl_ssu_ref_132.txt access2species_silva132.tsv > species2ncbiId.tsv.tmp

grep -v -e 'unculture\|metagenome\|Ambiguous_taxa' species2ncbiId.tsv.tmp > species_names_to_ncbi_id.tsv

more species_names_to_ncbi_id.tsv | sort | uniq | sort > species_names_to_ncbi_id.tsv.uniq

mv species_names_to_ncbi_id.tsv.uniq species_names_to_ncbi_id.tsv 

rm species2ncbiId.tsv.tmp species accessions cons_taxonomy_7_levels_silva_132_no_eukaryotes.txt

```

Using Silva v132 we ended up with 36,249 unique species names with their corresponding NCBI Taxonomy Ids. 

Ideally, *microbetag* would be able to include KO-annotated genomes for all this list of NCBI Taxonomy ids. 



