---
layout: default
title: Input files and settings
nav_order: 5
---

# Input files and settings
{: .no_toc }

---


## Input files

| File            | Description                                                    | requirement_status        |
|-----------------|----------------------------------------------------------------|---------------------------|
| abundance_table | An abundance table (in `.tsv`, `.csv`, `.biom` format)         | mandatory                 | 
| metadata_file   | File describing the sequencing data; for examples see [here]() | optionally with FlashWeave| 
| network_file    | if already built                                               | optional                  |


Please, make sure in case you provide your abundance table as a `.tsv` or `.csv` file that in the first column you have always the sequence identifier and a complete 7-level taxonomy in the last one.
*microbetag* requires for a 7-level taxonomy scheme; for example:

```bash
Bacteria;Firmicutes;Thermoanaerobacteria;Thermoanaerobacterales;Thermoanaerobacteraceae;Caldanaerobius;Caldanaerobius polysaccharolyticus
```

in case an entry reaches only to a higher taxonomic level, microbetag fills the entry with NA values

for example

```bash
Bacteria;Firmicutes;Thermoanaerobacteria;Thermoanaerobacterales;Thermoanaerobacteraceae
```

would become

```bash
Bacteria;Firmicutes;Thermoanaerobacteria;Thermoanaerobacterales;Thermoanaerobacteraceae;NA;NA;NA
```


{: .important-title}
> Curate your taxonomies! 
> 
> If you have a taxonomy that "skips" a level, or another one that has more levels, microbetag will fail. You need to curate those taxonomies manually and make sure you always have a 7-level scheme for all the entries on your table.




For input file examples, please have a look [here](https://github.com/hariszaf/microbetag/tree/develop/tests).


{: .warning } Do not use numeric characters only for labeling your samples and/or the sequences mentioned in your abundance table. For example, `324` as a sample id will lead microbetag to fail. 




### Case 1: all you have is your abundance table and your taxonomies 

For up to a few thousands of sequences entries, microbetag can be a one-stop-shop application performing both taxonomy annotation, network generation and annotation. 

Moving on with microbetag's taxonomy annotation is always a best practice as it makes sure that the sequences assigned to the species/strain level, they will all be annotated from microbetag. 

Further, one can keep both the taxonomies assigned from microbetag and from any other software. 

However, if you would like to move on with your taxonomy scheme, microbetag enables that but you should know that there's a big chance of loosing some annotations. 


In case 1, one may also have some metadata describing the sequencing data. FlashWeave, the software microbetag invokes to build the co-occurrence network, can exploit metadata. 

{: .important-title}
> METADATA FILE 
> 
> If you want to run FlashWeave with a metadata file, you need to remember that FlashWeave considers as variables both the sequence ids (i.e., ASVs/OTUs/bins) and the metavariables (e.g. pH, sex, any 
> variable on your metadata file). Thus, you need to have both of them as **rows**, contrary to what we do in most microbiome analyses. 
>
> Here is a toy example of how your files should look like: 
>
> `abundance_file.txt`
> 
> ```
> seqId    sample_1    sample_2    sample_3
> asv_1    10        0        3
> asv_2     0       21       43
> asv_3    32       31        2
> asv_4     0        0       12
>```
>
> `metadata_file.tsv`
>
>```
> Metadata_1      0.2     1.7       0
> Metadata_2      Yes      No       Yes 
>```
> As shown, the sample names are omitted from the `metadata_file.tsv`. 
> You need to make sure that their corresponding values are in the exact same order as in the `abundance_file.txt`. 
> In case the files are not provided like this, microbetag and/or the Docker image of microbetag preprocess, will fail.


{: .important-title}
> ADVANCED USAGE
>
>If you would like to have extra arguments for FlashWeave, then all you need to do is to run the `prep` image interactively and edit the `flashweave.jl` script accordingly (see [below](#the-preparation)). 




### Case 2: you already have a co-occurrence network 

In this case, you need to provide microbetag with both your abundance table and the co-occurrence network file. 
The latter can be of any form if microbetag is performed through the CytoscapeApp.
Otherwise, please make sure you provide the network as an edge file, for instance: 

{: .note}
>ASV_963239	ASV_4372091	0.3769868016242981
>
>ASV_4480529	ASV_4472202	0.4468387961387634
>
>ASV_4472202	ASV_4374302	0.4154910147190094
>
>ASV_4480529	ASV_4439469	0.39721810817718506


### Case 3: you like things the microbetag way

To get the optimal annotations in the more robust way, we **strongly suggest** you first prepare your data using the `prep_microbetag` Docker/Singularity image.
That will be almost always the case when you have large datasets with more than a few thousands of sequences and no network for them. 
Yet, even if you have a network, it is still **strongly suggested** to run the *taxonomy assignment* step, so microbetag can map more efficiently the taxa present to their corresponding GTDB genomes. 

Have a look at the ["preparation"](#the-preparation) section for how to do so! 



## Settings


There are cases where GTDB or Silva might have been used but with different formats than their original. 
For example, DADA2 has made a 


| Variable      | Description                       | Value |
|---------------|-----------------------------------|-------|
|`input_category`| If   | `abundance_table` \| `network` |
| `taxonomy` | In case a user's taxonomy is to be used, denotes which taxonomy scheme to be used from microbetag | [`GTDB` \| `dada2` \| `qiime2`]
| `phenDB`            | return phenotypic traits based on phen models  | bool |
| `faprotax`          | return annotations using the FAPROTAX database | bool |
| `pathway_complement`| return pathway complmementarities between associated nodes | bool |
| `seed_scores`       | return complementarity and cooperation scores based on metabolic reconstructions seed sets | bool |
| `manta`             | return clusters of nodes on the network using the manta package | bool | 
| `get_children`      | use genomes of children taxa of the taxa in the abundance table based on the NCBI Taxonomy scheme | bool |
| `heterogeneous`     | (FlashWeave) enable heterogeneous mode for multi-habitat or -protocol data with at least thousands of samples (FlashWeaveHE)| bool | 
| `sensitive`     | (FlashWeave) enable fine-grained associations (FlashWeave-S, FlashWeaveHE-S), sensitive=false results in the fast modes FlashWeave-F or FlashWeaveHE-F | bool | 



## The "preparation" 


Microbetag is a one-stop-shop application as it supports the taxonomical annotation of ASVs/OTUs, the building of the co-occurrence network and 
its annotation. 
However, its main goal is the latter and at the same point, the first two tasks can be computationally expensive especially for large datasets. 
To this end, a Docker/Singularity image is available supporting the taxonomy assignment of the ASVs/OTUs with GTDB taxonomies 
[a taxonomy annotated abundance file with the 16S GTDB (v.207) taxonomies](https://zenodo.org/records/6655692) and the creation of the co-occurrence network if asked. 


[Docker](https://docs.docker.com/get-docker/) or [Singularity](https://docs.sylabs.io/guides/3.0/user-guide/installation.html) needs to be installed. 
Then donwload the `prep_microbetag` image either by running: 


```bash
docker pull hariszaf/prep_microbetag
```


or 
```bash
singularity pull docker://hariszaf/prep_microbetag
```


Then you need to also download the<a href="https://github.com/hariszaf/microbetag/raw/preprocess/preprocess/test/config.yml" download="config.yml">`config` file</a> and edit it accordingly. 


{: .important}
The `cofing.yml` file and the input files to be used, need to be in the directory to be mounted (see Docker and Singularity commands below).


### I/O folder 

The folder to be mounted needs two mandatory files:
- an abundance table 
- the `config.yml` file


Two main steps are to be performed from this preparation image. 
A taxonomy annotation of the OTUs/ASVs to GTDB taxonomies using the IDTAXA algorithm of the DECIPHER package and the 16S sequences of the GTDB genomes or/and the building of a co-occurrence network using FlashWeave. 

The user can select which tasks to run through the `16s_gtdb_taxonomy_assign` and the `build_network` arguments. 

A thorough description of each argument can be found below as well as in the `config.yml` template.


|**Parameter**                |**Description**                                                                                         |
|-----------------------------|--------------------------------------------------------------------------------------------------------|
|``abundance_table_file``     | An OTU/ASV abundance table with a sequence identifier in first column and the sequence in the last one |
|``metadata_file``            | A metadata file with samples as rows and descriptions characterizing the samples as columns For an example file see https://github.com/meringlab/FlashWeave.jl/tree/master/test/data/HMP_SRA_gut. |
|``build_network``            | Build co-occurrence network using FlashWeave |
|``16s_gtdb_taxonomy_assign`` | Taxonomy assignment of the sequences using IDTAXA and the 16S sequences of the GTDB genomes. |
|``flashweave_sensitive``     | (From FlashWeave documentation) Enable fine-grained associations (FlashWeave-S, FlashWeaveHE-S), sensitive=false results in the fast modes FlashWeave-F or FlashWeaveHE-F  |
|``flashweave_heterogeneous`` | (From FlashWeave documentation) Enable heterogeneous mode for multi-habitat or -protocol data with at least thousands of samples (FlashWeaveHE) |
|``output_directory``         | Output folder name; it will be created within the mounted folder |


Now, based on your container technology you are ready to run the preparation image.

An example of a directory to mount can be seen [here](https://github.com/hariszaf/microbetag/tree/preprocess/preprocess/test). 
The mandatory abundance table file can be provided as a `.tsv` or a `.csv` file and needs to be specified in the `config.yml` file accordingly.


### Docker

To run directly 

```bash
docker run --rm -v ./test/:/media hariszaf/prep_microbetag
```

in this case it is the `./test/` local directory you mount in the `/media` folder on the container.


If you would like to initiate an interactive container you may run: 

```bash
docker run --rm -it --entrypoint /bin/bash -v ./test/:/media prep_microbetag
```
this would initiate a console from within the container with you as a root user:

```bash
root@69bdeedb582b:/pre_microbetag#
```


{: .highlight}
This case can be useful when several FlashWeave arguments not included in the basic config file need to be edited. 

You can see what is present under the `/pre_microbetag` directory:
```bash
root@69bdeedb582b:/pre_microbetag# ls
classify.R  flashweave.jl  gtdb_16s.RData  prep.py
```

and edit scripts, e.g.:

```bash
vim flashweave.jl
```

### Singulariity

The equivalent commands in Singularity would be :

To run the preprocess directly:
<!-- --env USER_ID=$HOST_USER_ID --env HOST_GID=$HOST_GROUP_ID -->
```bash
singularity run -B ~/prep_test/:/media prep_microbetag_latest.sif 
```

where again, you mount (`-B`) the `~/prep_test/` local directory to the `/media` folder on the container.

Likewise, to open a console:

```bash
singularity exec -B ~/prep_test/:/media prep_microbetag_latest.sif bash
cd /pre_microbetag/
```


