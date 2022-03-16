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


