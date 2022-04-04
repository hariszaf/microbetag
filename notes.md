# Sketching `microbetag`


## The workflow

- [FlashWeave](https://github.com/meringlab/FlashWeave.jl) [1] (if user does not provide an edge list along with the OTU table)

- [FAPROTAX](http://www.loucalab.com/archive/FAPROTAX/lib/php/index.php) (as part of the phenotipic channel)


> *microbetag* needs to assign a NCBI Taxonomy Id to those taxonomies of the OTU table that are at the species level. This is not a straight-forward task across the various reference databases might be used for the taxonomy assignment step. Currently, *microbetag* supports OTU tables derived from the Silva v.138 release; Qiime2 and DADA2 . Alternatively, one may provide an OTU table derived from any reference database of choice with an extra column denoting the NCBI Taxonomy Id assigned to each OTU (see example [here]()). 


### I/O 

**Mandatory** input files:

- OTU table 
- `config.yml` file

**Optional** input files:

- co-occurrence network 

- metadata table


**Output** files:

- a network file (formats to be supported: `.json`, `.gml`, `.csv`) file with the annotations from all the different modules in each node (taxon) and edge (association) as well as a notion of whether the association can be considered as *of higher confidence*.


## Approach 

It is qute often that an OTU table has many non-species level taxonomies. 
As a result, in microbial co-occurrence networks, there are associations 
both at the species but also in higher levels. 

**microbetag** intends to address this challenge by exploiting different sources of information depending on the taxonomic level 
of the taxa involved in each association pair. 

<!-- > Future work: BugBase is supposed to work with shotgun metagenomics as well. 
If that is so, then it would be great for *microbetag* to go that way. 
In any case, annotating associations coming from microbial co-occurrence networks derived from shotgun metagenomics can be of higher confidence in general; 
especially when they correspond to MAGs  -->




1. amplicon data: 

   * if a taxonomy does not exist, the OTUs/ASVs sequences are required. *microbetag* will use the phylogeny tree built out of the 16S sequences derived from the GTDB genomes (files: `ar122_ssu_reps_r202.fna` and `bac120_ssu_reps_r202.fna` files from the [GTDB FTP](https://data.gtdb.ecogenomic.org/releases/release202/202.0/genomic_files_reps/)) to place user's OTUs/ASVs on it before the annotation step (suggested in terms of the quality of the results)

   * if a taxonomy exists, then this needs to be based on Silva and with tab-separated columns, one column per sample and one row per OTU, using the 7-level format.


2. shotgun data: 

   * if a taxonomy does not exist the user needs to provide the genomes. *microbetag* will then run GTDB-tk to assign them, this will take some time. Alternatvely, the user may provide only the 16S sequences of his/her genomes and *microbetag will do as in the amplicon case.

   * if a taxonomy table is provided, it needs to have come using GTDB and the accession number of the closest genome to be provided too (suggested in terms of computation time)


In any case, *microbetag* will use the NCBI Taxonomy Id assigned to the species and strains present to extract relative infromation from the phenotypical and environmental data module, and a genome accession id to use the closest GTDB genome to each node for the pathway complementarity module. 




## Module 1 : Pathway complementarity

A strategic point of contradiction in the pathway complementarity module, is that 
*microbetag* focuses mostly in microbial associations networks coming from 16S rRNA amplicon data. 
Thus, it is most probably to get taxa at the genus level, or even at the species level. 
However, genomes are at the strain level. 
So we need to find a way to deal with this. 

At the same time, NCBI Taxonomy has been commonly used for a long time. 
However, GTDB Taxonomy is gaining space regarding the genome and the MAGs submission. 
So, the choice of the taxonomy about to be used is another strategic decision. 


### KEGG

<!-- KEGG GENOMES are **manually curated** and this makes their contribution essential to our task. 

- From KEGG organisms 
   - Bulk download - needs to pay first
   - API, for example: http://rest.kegg.jp/link/ko/miu

   The organisms in this case is denoted as "miu" that stands for Mitsuaria
   One can **map** the NCBI Taxonomy Ids to KEGG Ids by making use of this [file](https://www.genome.jp/kegg-bin/download_htext?htext=br08610&format=htext&filedir=); for browsing on the file click [here](https://www.kegg.jp/brite/br08610) instead.
   In the last line of this file you may see when it was last updated 
   ```
   #Last updated: Feburary 15, 2022
   ```
   To get `KEGG_IDs : NCBI Taxonomy` IDs pairs, you just need to run:

   ```
   grep -A 1 "\[TAX:" br08610.keg > pairs_in_2_lines
   awk -F" " '{print $2}' pairs_in_2_lines > KEGG_IDS
   sed -i '/^[[:space:]]*$/d' KEGG_IDS
   grep "\[TAX:" br08610.keg  > NCBI_IDS
   sed -i 's/.*\[TAX://g ; s/\]//g'  NCBI_IDS 

   ```

> **WRONG!**
The way described above is not correct! 
This is why there are more than 1 KEGG ids that are linked to the same 
NCBI Taxonomy Id!  -->



**Logical Expressions**

The pathway module is defined by the logical expression of K numbers, and the signature module is defined by the logical expression of K numbers and M numbers, allowing automatic evaluation of whether the gene set is complete, i.e., the functional unit is present, in a given genome. A space or a plus sign, representing a connection in the pathway or the molecular complex, is treated as an AND operator and a comma, used for alternatives, is treated as an OR operator. A minus sign designates an optional item in the complex.

Each space-separated unit is called a block, and the distinction is made for:
- complete modules
- incomplete but almost complete modules with only 1 or 2 blocks missing
- all modules that contain any matching K numbers

when evaluating the completeness check, such as in [KEGG Mapper](https://www.genome.jp/kegg/mapper/).

The reaction module is defined by the logical expression of RC numbers (Reaction Class identifiers), but this expression is not currently used for any evaluation purpose.


It seems like when there is a "--" symbol in a module, there's a missing KO term for a reaction.


We assume that a module occurs in a species, if at least the 60% of the KO terms 
required are available in the species' genome. 

>**Remember!**
A certain KO term might contribute in more than 1 KEGG modules ! 



### GTDB high quality representative genomes 

Using the latest version of GTDB ([v202](https://data.gtdb.ecogenomic.org/releases/release202/202.0/))
metadata file, we retrieved the ncbi genome accessions of the **representative** genomes 
of high quality, i.e. completeness > 95%  and contamination < 5%.

Using these accession numbers, we were able to download their corresponding `.faa` files when available. 

| genomes category | # of genomes|
|:-----:|:-------:|
|All GTDB  | 258,407 | 
| Representatives | 47,894 |
| High quality represntative | 26,778* | 
| HQ representative with `.faa` | 16,900** |

\* covering 22,009 unique NCBI Taxonomy Ids.

\** <REMEMBER TO COMPLETE THIS AS ABOVE>


We downloaded these 16,900 `.faa` files and then we used the 
[`kofamscan`](https://github.com/takaram/kofam_scan) software to annotate them with KO terms; `kofamscan` is a KEGG-family software. 
For more about it, you may see at the [KEGG website](https://www.genome.jp/tools/kofamkoala/) or/and the corresponding [publication](https://doi.org/10.1093/bioinformatics/btz859).

> Between annotations coming from the KEGG genomes and those from MGnify catalogues and GTDB high quality representative genomes there is a crucial difference. The first are manually curated, meaning their annotations have a greater level of confidence. *microbetag* combines these sets to cover the more taxa with the highest level of confidence possible, indicating which terms are coming from which genomes. 



* ssu_all_<release>.tar.gz
        FASTA file containing 16S rRNA sequences identified across the set of GTDB genomes passing QC. The assigned taxonomy
        reflects the GTDB classification of the genome. Sequences are identified using nhmmer with the 16S rRNA model (RF00177)
        from the RFAM database. All sequences with a length >=200 bp and an E-value <= 1e-06 are reported.
* bac120_ssu_reps_<release>.tar.gz
        FASTA file of 16S rRNA gene sequences identified within the set of bacterial representative genomes. The longest
        identified 16S rRNA sequence is selected for each representative genomes. The assigned taxonomy reflects the
        GTDB classification of the genome. Sequences are identified using nhmmer with the 16S rRNA model (RF00177) from the
        RFAM database. Only sequences with a length >=200 bp and an E-value <= 1e-06 are reported. In a small number of cases,
        the 16S rRNA sequences are incongruent with this taxonomic assignment as a result of contaminating 16S rRNA sequences.





## Module 2: Phenotypical data

### FAPROTAX

From the FAPROTAX website:
*The FAPROTAX database is optimized for taxonomies in SILVA releases 128 and 132. Other taxonomies not consistent with SILVA 132 may work sub-optimally.*


In the *microbetag* framework, the taxonomies assigned are first filtered and only those from the family level and lower are kept. 
Then FAPROTAX will run for the rest of the taxonomies of the OTU table. 



We had to edit the script.
Line 208 was changed to 
```python=
return (s.lower() != 'nan' and is_number(s))
```

```bash=
# working
./collapse_table.py -i finalTable.tsv -o func_table.tsv -g FAPROTAX.txt -d "Classification" -c "#" -v

# or to get OTU table per functional category 

./collapse_table.py -i finalTable.tsv -o func_table.tsv -g FAPROTAX.txt -d "Classification" -c "#" -v --out_sub_tables_dir .

```

### Traits from all over 

Super useful [repo](https://github.com/bacteria-archaea-traits/bacteria-archaea-traits). 

Relative publication at [Scientific data](https://www.nature.com/articles/s41597-020-0497-4) journal. 

The repo exploits both the NCBI and the GTDB taxonomies.

However, it is the NCBI Taxonomy ids that are used as keys. 
Moreover, a certain NCBI Taxonomy id may occur more than one in the `condensed` files 
and that is the `traits` files but only once in the `species` files. 

*microbetag* will use the `species` file if it is a species under study 
and the `traits` if it is a strain. 
If a strain is not present in the `traits` file, then the species will be used. 




### PhenDB 

https://phendb.org/




All the RefSeq genomes will be soon available for *microbetag* using the by default scores.



## Resources & approaches we `microbetag` could integrate 

Future work as in bullet-points:

- Exception case: time series data
- Exploit [EnDED](https://github.com/InaMariaDeutschmann/EnDED) tool
- Exploit [Metage2Metabo](https://github.com/AuReMe/metage2metabo) software
- The [equilibrator](https://equilibrator.weizmann.ac.il/) library 
- GUI



## References 


[1] Tackmann, Janko, João Frederico Matias Rodrigues, and Christian von Mering. "Rapid inference of direct interactions in large-scale ecological networks from heterogeneous microbial sequencing data." Cell systems 9.3 (2019): 286-296.

[2] Parks, Donovan H., et al. "A complete domain-to-species taxonomy for Bacteria and Archaea." Nature biotechnology 38.9 (2020): 1079-1086.

[3] Kanehisa, Minoru, et al. "KEGG: new perspectives on genomes, pathways, diseases and drugs." Nucleic acids research 45.D1 (2017): D353-D361.

[4] Mitchell, Alex L., et al. "MGnify: the microbiome analysis resource in 2020." Nucleic acids research 48.D1 (2020): D570-D578.

[5] Louca, Stilianos, Laura Wegener Parfrey, and Michael Doebeli. "Decoupling function and taxonomy in the global ocean microbiome." Science 353.6305 (2016): 1272-1277.

[6] Ward, Tonya, et al. "BugBase predicts organism-level microbiome phenotypes." BioRxiv (2017): 133462.



<!-- 
## Visualization libraries and alternatives 

### What about DASH Cytoscape? 

#### Dash Cytoscape as Docker app 


We followed  a [`plotly` community thread](https://community.plotly.com/t/running-dash-app-in-docker-container/16067) to build a web app of our own. 

An example of a Docker app is described [`here`](https://docs.docker.com/compose/gettingstarted/) at the Docker communication. 

We need to consider and study about: 

- [`Flask`](https://palletsprojects.com/p/flask/), a light Web Server Gateway Interface web application framework designed to make getting started quick and easy, with the ability to scale up to complex applications. 

- [`redis`](https://www.fullstackpython.com/redis.html), an in-memory key-value pair database typically classified as a NoSQL database. Redis is commonly used for caching, transient data storage and as a holding area for data during analysis in Python applications.


 -->









<!-- 
### Learn `neo4j`, `Cypher` and more (probably not to use at the time)

### Cypher Query Language


Cypher is a graph query language that is used to query the Neo4j Database. 
Cypher is like SQL - a declarative, textual query language.

In writing a cypher query, a **node** is enclosed between a **parenthesis** — like `(p:Person)` where `p` is a variable and `Person` is the type of node it is referring to.

**Relationships** are enclosed in square **brackets** - like `[w:WORKS_FOR]` where `w` is a variable and `WORKS_FOR` is the type of relationship it is referring to.

Putting it all together, the query below will return all people who have worked in Cloud Atlas movie.

```cypher=
MATCH (p:Person)-[relatedTo]-(m:Movie {title: "Cloud Atlas"})
RETURN p, m, relatedTo
```

Cypher is easy to learn and even more powerful and concise to write than SQL.

In graph, because we don't need join tables, Cypher expresses connections as graph patterns and is easy to read.

For example, if you want to find all actors that Tom Hanks has worked with, then below are the Cypher and SQL queries to fetch the data.

Cypher
```cypher=
MATCH (tom:Person {name: 'Tom Hanks'})-[a:ACTED_IN]->(m:Movie)<-[rel:ACTED_IN]-(p:Person)
return p, a, rel, m, tom
```

SQL
```sql=
SELECT * FROM persons p JOIN roles r ON (p.id = r.person_id)
JOIN movies m ON (m.id = r.movie_id)
JOIN roles r2 ON (m.id = r2.movie_id)
JOIN persons p2 ON (p2.id = r2.person_id)
WHERE p2.name = "Tom Hanks"
```

For more about Cypher, you may have a look [here](https://neo4j.com/developer/cypher/).


### AuraDB

*What are the the AuraDB Free database limits?*

The free database is limited on graph size: 50K nodes and 175K relationships.


*How long can I use the Free database?*
The database is free forever, no credit card required. 
It will **be automatically paused** if you do not use it for **72 hours**. But don't worry, you can resume the database at any time within 90 days.


### neo4j 

neo4J is a **native** graph database. 
Each node contains direct pointers to all the nodes that it is connected to.
These direct pointers are called *relationships*.
All the information needed to find the next node in a sequence is available in the node itself. 
The native storage layer is a connected graph, that's what *native* means.

Because of this principle neo4j does not need to compute the relationships between your data at query time, 
the connections are already there stored right in the database. 
Therefore, neo4J does not need an index system!




 -->
