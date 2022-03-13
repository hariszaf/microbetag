# Sketching `microbetag` - building new skills!


## `microbetag` scheme

> Taxonomies supported : NCBI Taxonomy Vs GTDB Taxonomy


**Input points:**

- OTU table & metadata table

- Co-occurrence network & metadata table

- Exception case: time series data

**Pipeline*: **

- FlashWeave (if user does not have a co-occurrence network as input but an OTU table)

- EnDED (optional; if yes run 1st)

- BUGBASE \& FAPROTAX (for the phenotipic channel)

- Pathway complementarity module 

- EcoEFMSs

\* I am thinking that as all these modules are actually different ways
to approach our goal, we could ask the user to build his/her own pipeline 
and not having a default one! 
Having a default makes things easier sure, but maybe in this case it would be better
to ask the user to *build* the pipeline. 

**Output:** 

- `.tsv` files with the annotations per edge 

- a GUI with the annotated graph



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

KEGG GENOMES are **manually curated** and this makes their contribution essential to our task. 

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
NCBI Taxonomy Id! 

>**Remember!**
A certain KO term might contribute in more than 1 KEGG modules ! 


**Logical Expression**

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




### MGnify catalogues





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


## Module 2: Phenotypical data

### FAPROTAX

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








## Resources & approaches we `microbetag` could integrate (future work)

- https://equilibrator.weizmann.ac.il/





## 






## Visualization libraries and alternatives 

### What about DASH Cytoscape? 

#### Dash Cytoscape as Docker app 


We followed  a [`plotly` community thread](https://community.plotly.com/t/running-dash-app-in-docker-container/16067) to build a web app of our own. 

An example of a Docker app is described [`here`](https://docs.docker.com/compose/gettingstarted/) at the Docker communication. 

We need to consider and study about: 

- [`Flask`](https://palletsprojects.com/p/flask/), a light Web Server Gateway Interface web application framework designed to make getting started quick and easy, with the ability to scale up to complex applications. 

- [`redis`](https://www.fullstackpython.com/redis.html), an in-memory key-value pair database typically classified as a NoSQL database. Redis is commonly used for caching, transient data storage and as a holding area for data during analysis in Python applications.












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
