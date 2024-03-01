# microbetag

`microbetag` attempts to be a microbial interactions co-occurrence network annotator

In this branch we have a stand-alone tool to run the software on your own bins.
This module cannot be performed directly from the CytoscapeApp as it requires intense computing resources is not possible to provide in a web-based app.
Yet, once you run this stand-alone version, you will still be able to load the annotated network returned on Cytoscape and investigate the annotations through the MGG CytoscapeApp of ours.

## Dependencies

- Docker
- kofam_database



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
python microbetag.py -conf config.yml 
```

> **Taxonomy**
>
> Contrary to the microbetagDB - dependent approach, in this case we do not care on the scheme!
>


## Docker

The `microbetag` annotator will be available as a Docker image....

## Graphical User Interphase (GUI)

But it `microbetag` will come with a GUI too.

A first thought on how to do this is to follow this scheme:

- build a database with the `microbetag` annotations
- implement `microbetag` so it asks queries on the db
- output could be a list of files, some of them `.html` so they could provide interactive visualizations such as having annotations on the nodes or/and edges

## Funding

This project is funded by an [EMBO Short-Term Fellowship](https://www.embo.org/funding/fellowships-grants-and-career-support/scientific-exchange-grants/).
