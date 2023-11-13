---
layout: default
title: (optimal) input files and settings
nav_order: 5
---

# Input files and settings
{: .no_toc }

---


## Input files

In case a co-occurrence network is already available, 
microbetag wi annotate it. 

If there is 




## Settings


Input category: `[abundance table | network]`


Taxonomy: `[GTDB | Silva | other]`

There are cases where GTDB or Silva might have been used but with different formats than their original. 
For example, DADA2 has made a 




## The "preparation" 


Microbetag is a one-stop-shop application as it supports the taxonomical annotation of ASVs/OTUs, the building of the co-occurrence network and 
its annotation. 
However, its main goal is the latter and at the same point, the first two tasks can be computationally expensive especially for large datasets. 
To this end, a Docker/Singularity image is available supporting the taxonomy assignment of the ASVs/OTUs with GTDB taxonomies 
[a taxonomy annotated abundance file with the 16S GTDB (v.207) taxonomies](https://zenodo.org/records/6655692) and the creation of the co-occurrence network if asked. 


[Docker](https://docs.docker.com/get-docker/) or [Singulariity](https://docs.sylabs.io/guides/3.0/user-guide/installation.html) needs to be installed. 
Then donwload the `prep_microbetag` image either by running: 


```bash
docker pull hariszaf/prep_microbetag
```


or 
```bash
singularity pull docker://hariszaf/prep_microbetag
```


Then you need to also download the<a href="https://github.com/hariszaf/microbetag/raw/preprocess/preprocess/test/config.yml" download="config.yml">`config` file</a> and edit it accordingly. 


To have ownership and permissions to use the data products of the container, please first keep the following two variables:

```bash
HOST_USER_ID=$(id -u)
HOST_GROUP_ID=$(id -g)
```

and then based on your container technology you are ready to run the preparation image.

### Docker

To run directly 

```bash
docker run --rm -v ./test/:/media -e USER_ID=$HOST_USER_ID  -e HOST_GID=$HOST_GROUP_ID   hariszaf/prep_microbetag
```

If you would like to initiate an interactive container you may run: 

```bash
docker run --rm -it --entrypoint /bin/bash -v ./test/:/media -e USER_ID=$HOST_USER_ID  -e HOST_GID=$HOST_GROUP_ID   prep_microbetag
```

{: .highlight}
This case can be useful when several FlashWeave arguments not included in the basic config file need to be edited. 


### Singulariity

The equivalent commands in Singularity would be :

To run the preprocess directly:
```bash
singularity run -B ~/prep_test/:/media --env USER_ID=$HOST_USER_ID --env HOST_GID=$HOST_GROUP_ID prep_microbetag_latest.sif 
```

and likewise, to open a console:

```bash
singularity exec -B ~/prep_test/:/media --env USER_ID=$HOST_USER_ID --env HOST_GID=$HOST_GROUP_ID prep_microbetag_latest.sif bash

cd /pre_microbetag/

```


### `config.yml` file for the preparation


Two main steps are to be performed from this preparation image. 
A taxonomy annotation of the OTUs/ASVs to GTDB taxonomies using the IDTAXA algorithm of the DECIPHER package and the 16S sequences of the GTDB genomes or/and the building of a co-occurrence network using FlashWeave. 

The user can select which tasks to run through the `16s_gtdb_taxonomy_assign` and the `build_network` arguments. 

A thorough description of each argument can be found below as well as in the `config.yml` template.


|**Parameter**                |**Description**                                                                                         |
|-----------------------------|--------------------------------------------------------------------------------------------------------|
|``abundance_table_file``     | An OTU/ASV abundance table with a sequence identifier in first column and the sequence in the last one |
|``metadata_file``            |  Using the filtered and merged sequences, it returns a taxonomic inventory                             |
|``build_network``            |  Exports coding sequences                                                                              |
|``flashweave_sensitive``     |  Performs functional annotation on the coding genes found using a list of resources: InterPro, KEGG    |
|``flashweave_heterogeneous`` |  Assembles the filtered and merged sequences to contigs                                                |
|``output_directory``         |  Assembles the filtered and merged sequences to contigs                                                |
|``16s_gtdb_taxonomy_assign`` |  Assembles the filtered and merged sequences to contigs                                                |

