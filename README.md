# microbetag

`microbetag` attempts to be a microbial interactions co-occurrence network annotator

In this repository we have:

- the precalculation steps/scripts under the `microbetagDB` folder; the database of `microbetag` consists of the outcome of these precalculations and to get those the resources under the `ref-dbs` folder are used
- the `microbetag` python library which makes use of the microbetagDB and assumes several other software have been installed properly (see Installation)
- a `Dockerfile` to run `microbetag` in a container-based way

## Dependencies

`microbetag` makes use of several 3rd-party software. Thus, it inherits the dependencies of theirs.

* Julia
* FlashWeave
* FAPROTAX
* BugBase

To help you create such an environment, we provide the [`installation.sh`](installation.sh) script. However, this will work smoothly in case you are running `microbetag` on Linux systems. If your OS is MacOs or Windows, you need to set this manually or jump into the containerized version of `microbetag`.

<!-- 

MUTED

BugBase requires the `biomformat` library. To install it on R version 4.2 you may run:

```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("biomformat")

```

-->


## Installation

## How to run

```bash
python microbetag.py -conf config.yml 
```

`microbetag` expects a **7-level taxonomy** in a `.csv` or `.tsv` format with:

* **no** `D_1__` or `K__`,`P__` any other similar prefixes befor the taxon name
* the complete species name if that is available on its last column.

> **Example**
>
> If your taxonomy is
>
> `Bacteria;Firmicutes;Thermoanaerobacteria;Thermoanaerobacterales;Family III;Thermoanaerobacterium;thermosaccharolyticum`,
>
> you need to set it as:
>
> `Bacteria;Firmicutes;Thermoanaerobacteria;Thermoanaerobacterales;Family III;Thermoanaerobacterium;Thermoanaerobacterium thermosaccharolyticum`

One

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
