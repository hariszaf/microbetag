# Mapping files for the on-the-fly version

microbetag's on-the-fly version, needs two types of mapping files compared to the stand-alone case:

	- mapping taxonomies to their NCBI Taxonomy IDs, those to their corresponding GTDB representative genomes, and the last ones 
	to their corresponding PATRIC ids, available on microbetagDB

	- for the seed complementarity step, microbetag needs to make use of the pre-calculated seed and and non-seed sets of the PATRIC genomes available. 
	Also, since we focus on KEGG MODULES, we need to match which of those compounds do participate in KEGG MODULES.

In this folder, we keep mapping files for the on-the-fly version only.

Here comes their short description:



|                File                 | Description  |
|:-----------------------------------:|:------------:|
|`kegg_maps_descriptions.tsv`         |A 2-column tab separated file with KEGG PATHWAYS (maps) and their description | 
|`module_descriptions`                |A 3-column tab separated file with the metabolism category of a KEGG MODULE, the module ID and its description per line. | 
|`related_kegg_maps_descriptions.tsv` | A 3-column tab separated file with the KEGG PATHWAY (map) ID, its name and the KEGG metabolism category it is part of| 
|`seedId_keggId_module.tsv`           | A 3-column tab separated file with the ModelSEED compound IDs, their corresponding one in KEGG and a KEGG MODULE it is part of. Attention! A compound may participate in more than one module. In this case, several entries for the compound will be recorded.| 



And for taxonomy: 

|                File                   | Description  |
|:-------------------------------------:|:------------:|
|`species2ncbiId.tsv`                   | A 2-column tab separated file with NCBI species name with their corresponding NCBI Taxonomy Ids |
|`genera2ncbiId.tsv `                   | A 2-column tab separated file with the NCBI genus name and its corresponding NCBI Taxonomy ID   |
|`families2ncbiId.tsv`                  | A 2-column tab separated file with the NCBI family name and its corresponding NCBI Taxonomy ID  |
|`species2ncbiId2accession.tsv`         | A 3-column tab separated file with species name, their corresponding NCBI Taxonomy IDs and their assigned GTDB representative genomes. Attention! A species name may map to several NCBI Taxonomy IDs and apparently to several genomes. To be used when taxonomy scheme is `Ohter` in the on-the-fly version. |
|`gc_accession_16s_gtdb_ncbid.tsv`      | A 3-column tab separated file with species names, NCBI Taxonomy IDs and their corresponding genomes for the case where taxonomies were assigned using the `microbetag_prep` tool. In this case, a GTDB v207, DADA2 formatted 16S rRNA gene sequences for both bacteria & archaea is being used https://zenodo.org/records/6655692. |
|`gtdbSpecies2ncbiId2accession.tsv`     | Like the above, but with the original GTDB species names. |
|`gtdb_silvaSpecies2ncbi2accession.tsv` | Like the above, but using the Silva species names. |


The following files, should be now deprecated:

- `dada2ncbi2accession.tsv`
- `qiime2species2ncbi2accession.tsv`

