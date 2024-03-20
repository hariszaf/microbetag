# microbetag

`microbetag` attempts to be a microbial interactions co-occurrence network annotator

In this branch we have a stand-alone tool to run the software on your own bins.
This module cannot be performed directly from the CytoscapeApp as it requires intense computing resources is not possible to provide in a web-based app.
Yet, once you run this stand-alone version, you will still be able to load the annotated network returned on Cytoscape and investigate the annotations through the MGG CytoscapeApp of ours.

## Dependencies

- Docker
- kofam_database
- gurobi license if you want to use carveme 
    https://support.gurobi.com/hc/en-us/community/posts/4406485885841-Installing-Gurobi-on-a-Docker-container-Ubuntu

Attention! You will need a Web License Manager not a typical license.

```bash
mkdir kofam_database &&\
    cd kofam_database &&\
    wget -c ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz &&\
    wget -c ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz &&\
    gzip -d ko_list.gz &&\
    tar zxvf profiles.tar.gz 
```



## Installation

## How to run

```bash
docker run --rm -it  \
    --volume=./tests/dev_io_microbetag/:/data \
    --volume=./microbetagDB/ref-dbs/kofam_database/:/microbetag/microbetagDB/ref-dbs/kofam_database/ \
    --volume=$PWD/gurobi.lic:/opt/gurobi/gurobi.lic:ro \
    --entrypoint /bin/bash  microbetag

```

Issue where different findings from carveme using cplex and gurobi is mentioned:
https://github.com/cdanielmachado/carveme/issues/141#issuecomment-912309490


## Docker

The `microbetag` annotator will be available as a Docker image....

## Graphical User Interphase (GUI)

But it `microbetag` will come with a GUI too.

A first thought on how to do this is to follow this scheme:

- build a database with the `microbetag` annotations
- implement `microbetag` so it asks queries on the db
- output could be a list of files, some of them `.html` so they could provide interactive visualizations such as having annotations on the nodes or/and edges


<!-- 
modelseedpy:  
 "scikit-learn == 1.2.0",   
 version lock for pickle ML models https://github.com/ModelSEED/ModelSEEDpy/blob/dev/setup.py
 "scipy >= 1.5.4",
 yet: python3 -m pip install scikit-learn==0.24.2 is the only working to load pickle.... 

 phenotrex: scikit-learn>=0.22,<=0.23.2   https://github.com/univieCUBE/phenotrex/blob/master/requirements/prod.txt
-->


## Funding

This project is funded by an [EMBO Short-Term Fellowship](https://www.embo.org/funding/fellowships-grants-and-career-support/scientific-exchange-grants/).
