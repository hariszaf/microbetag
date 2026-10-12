# microbetag
_annotating networks with phenotypic traits and potential metabolic interactions_

## About

`microbetag` annotates microbial co-occurrence networks:
  - on the **node** level: with phenotypic traits of the representing taxa
  - on the **edge** level: with potential metabolic interactions among the co-occurrying or mutual exclusive taxa.


The simplest way to use `microbetag` is through its Cytoscape App, called `MGG`.
You can [download `MGG`](https://apps.cytoscape.org/apps/mgg) directly from the Cytoscape ApppStore.

`MGG` allows both to run `microbetag` on-the-fly and visualize the returned annotated network.
However, using the on-the-fly version, means `microbetag` will try to match your taxa to their closest GTDB representative genomes.
Then, it will annotate your network based on pre-calculations of those GTDB genomes available on `microbetagDB`.

To annotate your network using your own genomes, run microbetag locally.
The steps that run external tools on the genomes (gene prediction, KEGG annotation, trait prediction,
genome-scale model reconstruction, FAPROTAX, network inference and clustering) are computationally
expensive; they run in the [**microbetag-nf**](https://github.com/msysbio/microbetag-nf) Nextflow pipeline,
on your computer or an HPC. The `microbetag` package then computes the pathway and seed complementarities
from their outputs and builds the annotated network, e.g. `mtag_net_2025-05-08_17-47.cx2`.
Load this file in Cytoscape to go through its annotations as in the on-the-fly version.

| | Runs | Install |
|---|---|---|
| **microbetag-nf** | the external tools, then `microbetag` | Nextflow and Docker/Singularity; see [microbetag-nf](https://github.com/msysbio/microbetag-nf) |
| **microbetag** (this repository) | pathway and seed complementarities, the annotated network, the on-the-fly version | `pip install` |

If you already have the outputs of those tools, for example GEMs or a KEGG annotation of your genomes,
you can run `microbetag` on them directly.


For more on how to use `microbetag`, `MGG` and their features, feel free to have a look 
at [`microbetag`'s documentation page](https://microbetag.readthedocs.io/en/v1.0.3/).


## Installation: Running `microbetag`..

### Option 1: On-the-fly 

To use `microbetag` on-the-fly, all you have to do is to make sure you have Cytoscape on your computing system and add the `MGG` app.

* **You can [install Cytoscape](https://cytoscape.org/) from their website.**

* **After firing Cytoscape you can [add `MGG` by downloading it from the Cytoscape AppStore](https://apps.cytoscape.org/apps/mgg).**


### Option 2: The microbetag-nf pipeline

To go from your genomes to an annotated network, use the [microbetag-nf](https://github.com/msysbio/microbetag-nf)
pipeline. It runs the external tools in containers and then `microbetag`.


### Option 3: The microbetag package

`microbetag` needs Python 3.10 or later. From a clone of this repository:

```bash
git clone https://github.com/hariszaf/microbetag.git
cd microbetag
pip install .
```

Add `".[otf]"` for the on-the-fly version, which queries the microbetag database, and `".[test]"` to run the tests.

To check your installation, run the tests from the root of the repository:

```bash
pip install ".[test]"
pytest
```

[`tests/README.md`](./tests/README.md) describes how to run single tests, or single microbetag steps on the test data.

#### In a container

The `hariszaf/microbetag` image holds the package and its commands:

```bash
docker pull hariszaf/microbetag:<version>
```

or build it from the root of this repository:

```bash
docker build -t hariszaf/microbetag:<version> .
```

Mount the folder with your configuration YAML and the files it points to at `/data`.
Paths in the YAML are relative to the YAML's folder, so a configuration that works locally also works in the container:

```bash
docker run --rm -v $PWD/my_study:/data hariszaf/microbetag:<version> microbetag --config /data/config.yml
```

For Singularity / Apptainer, convert the image once (no `sudo` needed) and run it with the same mount:

```bash
singularity pull microbetag.sif docker://hariszaf/microbetag:<version>
singularity exec -B $PWD/my_study:/data microbetag.sif microbetag --config /data/config.yml
```


## Further `microbetag` links

Overall, `microbetag` is a software ecosystem composed of a collection of modular code components:

- `microbetag-nf`: the Nextflow pipeline running the external tools, and preparing amplicon data for the on-the-fly version; individual [GitHub repo](https://github.com/msysbio/microbetag-nf)
- `microbetagDB` : `db` branch of the [msysbio/microbetag GitHub repo](https://github.com/msysbio/microbetag/tree/db)
- `microbetag_prep`: replaced by the prep workflow of `microbetag-nf`; the former tool is on the `preprocess` branch of the [msysbio/microbetag GitHub repo](https://github.com/msysbio/microbetag/tree/preprocess)
- `microbetagApp`: individual [GitHub repo](https://github.com/msysbio/microbetagApp-public)
- `MGG` source code: individual [GitHub repo](https://github.com/hariszaf/MGG)
- pre-calculations: Seed and non-seed sets for GTDB representative genomes available in `microbetag`, their corresponding GEMs, models for `phenotrex`-based predictions and KEGG annotations of each of those genomes are available on [Zenodo](https://zenodo.org/records/15102937).
  
For potential issues when running `microbetag`, or ideas for further features, 
please feel free to join us on [`Matrix`](https://matrix.to/#/#microbetagcommunity:matrix.org).



## How to run `microbetag`

### On-the-fly

In case you are about to run `microbetag` using the on-the-fly option, you need to only provide you actual input data files, i.e.:
- an abundance file (mandatory)
- a co-occurrence network (optional)
- metadata file for inferring the co-occurrence network (optional)

For thorough descriptions on the format of each of these files, 
[continue on the `MGG` tutorials](https://microbetag.readthedocs.io/en/latest/basic_usage/mgg_totorials.html).


### Locally

`microbetag` needs a **configuration YAML file**. Each `microbetag` version comes with its own template,
which this command writes out:

```bash
microbetag --template config.yml
```

It points to:
- the **co-occurrence network**, as an edge list, and an **abundance table** or a sequence id to taxonomy map;
- for each annotation you want, the **output of the tool** it comes from:

| Annotation | Input | Produced by |
|---|---|---|
| Pathway complementarity | `ko_merged_file`: KEGG ORTHOLOGY terms per genome | kofam (hmmsearch) on the predicted genes |
| Seed complementarity | `genome_scale_models`: one GEM (.xml) per genome | CarveMe, gapseq or ModelSEEDpy |
| Phenotypic traits | `phenotrex_predictions` | phenotrex |
| FAPROTAX traits | `faprotax_sub_tables` | FAPROTAX |
| Network clusters | `prev_clustered_network` | manta |

microbetag-nf produces all of them. Switch off (`value: false`) the steps you have no input for.
Then run:

```bash
microbetag --config config.yml
```

[`tests/data/microbetag`](./tests/data/microbetag/) is a complete example, with its [`config.yml`](./tests/data/microbetag/config.yml).

The steps are also available as separate commands (`microbetag-pcompl`, `microbetag-seed-sets`,
`microbetag-annotate`, ...), which the microbetag-nf pipeline uses; each has `--help`.
`microbetag` can also be used as a Python library, e.g. to get the seed complementarities of a set of GEMs
without a network (see [`tests/unit/test_seed_compl.py`](./tests/unit/test_seed_compl.py)).


## RTD

`microbetag` uses ReadTheDocs and Sphinx for its documentation (see `docs/`).

> **For contributors**
>
> To test changes locally first, run:
>
> Make sure you install `requirements.txt` on a conda environment: 
> 
>```
> cd docs/
> pip install -r requirements.txt
>```
>
> and then, you may run:
> ```
> sphinx-build -b html -d _build/doctrees -D language=en . _build/html -v
> ```
> from within the `docs/` directory.
> 

A new RTD is being produced every time a new _microbetag_ `tag` is released automatically.



## Cite

If you use `microbetag`, please consider citing us:

Zafeiropoulos H, Michail Delopoulos EI, Erega A, Schneider A, Geirnaert A, Morris J, Faust K.
microbetag: simplifying microbial network interpretation through annotation, enrichment tests and metabolic complementarity analysis.
*Genome Biology* 26, 292 (2025). DOI: [10.1186/s13059-025-03769-2](https://doi.org/10.1186/s13059-025-03769-2)


## Funding

This project is funded by: 
- the [3D' omics](https://www.3domics.eu) Horizon project (101000309).
- an [EMBO Short-Term Fellowship](https://www.embo.org/funding/fellowships-grants-and-career-support/scientific-exchange-grants/)

