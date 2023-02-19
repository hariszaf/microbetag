# Linking GTDB to NCBI Taxonomy Ids

We are interested in mapping the species/strains present on the network into their corresponding NCBI Taxonomy Ids.


## Using the Qiime2 release of Silva 132 in 7-levels

Based on the `consensus_taxonomy_7_levels.txt` under the 99 folder of the taxonomy `.zip` file of the [Qiime files of Silva v.132](https://www.arb-silva.de/download/archive/qiime/), we build two-column files with the GTDB names of the 7-level taxonomy scheme and their corresponding NCBI Taxonomy ids. 


### Case 1


Get the `nodes.dmp` file of NCBI Taxonomy from its FTP ([taxdump.tar.gz](https://ftp.ncbi.nih.gov/pub/taxonomy/)). 

From this file we can export the NCBI Taxonomy Ids that correspond to each taxonomic level: 

```
grep "species" nodes.dmp | awk -F"\t" '{print $1}' > ncbi_species_ids.tsv

# likewise for the other 6 levels
```

then 

```bash
while IFS= read -r line; 
do 
awk -v event="$line" '$1==event {print $3"\t"$1}' names.dmp >> phyla.tsv ; 
done < ncbi_phyla_ids.tsv
```

Likewise, for all the taxonomic levels. 

This returns to something like the following: 

```bash
Bacteroides-Cytophaga-Flexibacter       976
Bacteroidetes   976
"Bacteroidetes" 976
Bacteroidota    976
Bacteroidota    976
"Bacteroidota"  976
BCF     976
CFB     976
CFB     976
```
meaning that for a NCBI Tax Id we have several names, synonyms etc. 

Then, we sort and keep unique entries by: 

```bash
sort classes.tsv | uniq | sort > classes.tsv.tmp
mv classes.tsv.tmp classes.tsv
```


Thus, we can now map 






### Case 2 

We use the `ncbi-taxonomist`

dada2 = pd.read_csv("silva_138_taxonomies.tsv", sep=";" , header = None, engine="python")
dada2.drop(7, axis=1, inplace=True)
dada2.to_csv("compl.tsv", sep="\t", header=False, index = False)

# Load names
ncbi_names_ids = pd.read_csv("/home/luna.kuleuven.be/u0156635/github_repos/microbetag/mappings/gtdb_ncbi/ncbi_dump/names.dmp", sep="\t|\t", header = None, engine = "python")


----------







Using the `taxonomies_per_tax_level.py` we come up with a one-column file per taxonomic level with the taxa names included. 

Then, we deal differently the species and the rest of the taxonomic levels.

For the species case, we run the 



### Phyla 

Remove `Ambiguous_taxa` 

There are cases with multiple hits. 
E.g. D_1__Aquificae  187857\n200783

Resolve manually. 


 
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


```
grep -v "Eukaryota" cons_taxonomy_7_levels_silva_132.txt > 
        qiime2_consensus_taxonomy_7_levels.tsv.txt

awk -F" " '{print $1}' qiime2_consensus_taxonomy_7_levels.tsv.txt > accessions

awk -F";" '{print $NF}' qiime2_consensus_taxonomy_7_levels.tsv.txt > species

paste -d "\t"  accessions species > access2species_silva132.tsv

./map_accession_to_ncbi_id.awk taxmap_embl_ssu_ref_132.txt access2species_silva132.tsv >  
        species2ncbiId.tsv.tmp

grep -v -e 'unculture\|uncultivated\|metagenome\|Ambiguous_taxa' species2ncbiId.tsv.tmp > 
        species_names_to_ncbi_id.tsv

more species_names_to_ncbi_id.tsv | sort | uniq | sort > species_names_to_ncbi_id.tsv.uniq

mv species_names_to_ncbi_id.tsv.uniq species_names_to_ncbi_id.tsv 

rm species2ncbiId.tsv.tmp species accessions qiime2_consensus_taxonomy_7_levels.tsv.txt

```

Using Silva v132 we ended up with 36,249 unique species names with their corresponding NCBI Taxonomy Ids. 

Ideally, *microbetag* would be able to include KO-annotated genomes for all this list of NCBI Taxonomy ids. 

In addition, running: 
```
awk -F";" '{print $5}' qiime2_consensus_taxonomy_7_levels.tsv.txt > family
```
we get all the entries at the family level. 

And, with 

```
awk -F";" '{print $6}' qiime2_consensus_taxonomy_7_levels.tsv.txt > genus
```

the genus. 


However, for those cases where GTDB uses names of lower taxonomic levels to the higher ones, we removed them from the latter ones. 
For example, *Virgulinella fragilis* is a eukaryota actually but its chloroplast is under the "Chloroplast" "family" of Silva. 
Its 7-level taxonomy on Silva would be: 
```
D_0__Bacteria;D_1__Cyanobacteria;D_2__Oxyphotobacteria;D_3__Chloroplast;D_4__Virgulinella fragilis;D_5__Virgulinella fragilis;D_6__Virgulinella fragilis
```
In this case, we skip of keeping that as a genus and a family name and we only keep it at the species level. 


To get the NCBI Taxonomy id for the non-species level taxa, we use the `ncbi-taxonomist` tool and the 




