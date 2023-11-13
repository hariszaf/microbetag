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


```bash
HOST_USER_ID=$(id -u)
HOST_GROUP_ID=$(id -g)
```



```bash
docker run --rm -it --entrypoint /bin/bash -v ./test/:/media -e USER_ID=$HOST_USER_ID  -e HOST_GID=$HOST_GROUP_ID   prep_microbetag
```

```bash
docker run --rm -it -v ./test/:/media -e USER_ID=$HOST_USER_ID  -e HOST_GID=$HOST_GROUP_ID prep_microbetag
```


singularity exec -B ~/prep_test/:/media --env USER_ID=$HOST_USER_ID --env HOST_GID=$HOST_GROUP_ID prep_microbetag_latest.sif /usr/bin/python3 /pre_microbetag/prep.py


