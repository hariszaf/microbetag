---
title: microbetag as a Nextflow pipeline
layout: default
parent: Additional tutorials
nav_order: 3
description: "running microbetag on your genomes with Nextflow"
---


# The pipeline: microbetag-nf

<div style="display: flex; gap: 10px;">
   <a href="https://github.com/hariszaf/microbetag/tree/develop/pipeline" class="btn-green" target="_blank"> pipeline </a>
</div>

The pipeline, **microbetag-nf**, takes you from your genomes to an annotated network.
It runs the external tools, each in its own container, and then the `microbetag` commands:

| Step | Tool | Gives `microbetag` |
|---|---|---|
| Gene prediction | Prodigal | |
| KEGG ORTHOLOGY annotation | kofam (hmmsearch) | the KO terms of each genome, for pathway complementarities |
| Genome-scale metabolic models | CarveMe, gapseq or ModelSEEDpy | the GEMs, for seed complementarities |
| Phenotypic traits | phenotrex | trait predictions |
| Literature-based traits | FAPROTAX | FAPROTAX sub-tables |
| Network inference, if you have no network | FlashWeave | the network |
| Network clustering | manta | the clustered network |

`microbetag` then computes the pathway and seed complementarities and builds the annotated network.
The pipeline spreads the work over the genomes, so it runs on a laptop as well as on an HPC.

The pipeline is one of the two ways to use `microbetag` (see [usage modes](../modes.md)), and it serves the other
one too:

- **with the package**, it annotates your networks with your own genomes: `precalc` once for the genomes, then
  `annotate` for each network; you browse the annotated network with `MGG`;
- **for a server**, the output of `precalc` becomes the [database of your own server](#your-own-database), which
  annotates networks sent from `MGG` with your genomes; the server itself runs the pipeline's `network` workflow
  for FlashWeave, FAPROTAX and manta;
- **for our server**, the `prep` workflow prepares large amplicon datasets: GTDB taxonomies for the sequences, and
  the network; see [the preparation steps](../tutorials_otf/prep.md).

If you already have the outputs of these tools, you can also run the package on them directly;
see the tutorial [.. using custom genomes](./local.md).

The pipeline is the `pipeline/` folder of the
<a href="https://github.com/hariszaf/microbetag" target="_blank">microbetag repository</a>;
paths in this page are relative to it.


## Overview with an example dataset

The `tests/data/genomes` folder of the repository contains 3 genomes in FASTA format; this is what we call a
_genome catalogue_.

* The `precalc.nf` workflow computes all pairwise relationships between the genomes of a catalogue.
* The `annotate.nf` workflow infers a network, if one is not provided, and uses these pre-calculations to annotate it.

The `tests/expected/precalc` folder has the pre-calculations obtained by running `precalc.nf` on the example genomes,
and `tests/expected/annotate` the results of `annotate.nf` with these pre-calculations and an example abundance
table. `tests/expected/annotate/test.cx2` is the annotated network; you may [load it in Cytoscape](../tutorials_core/load.md)
and go through it with the MGG app.

Take extra care with the [Nextflow configuration](#nextflow-configuration), since you may need to edit it for your
computing environment.


## Requirements

Workflows, sub-workflows and modules can all run with [Docker](https://www.docker.com/) or
[Singularity](https://sylabs.io/singularity/)/[Apptainer](https://apptainer.org/). You need:

1. [Nextflow](https://www.nextflow.io/docs/latest/install.html)
2. One of these containerization technologies, usually depending on the computing system:
   - local PC/laptop: [Docker](https://www.docker.com/get-started/)
   - HPC: [Singularity](https://docs.sylabs.io/guides/3.0/user-guide/installation.html)/[Apptainer](https://apptainer.org/docs/admin/main/installation.html)

For [pathway complementarity](../modules/modules.md), if you do not already have a KEGG
ORTHOLOGY annotation of your genomes, you need the KOfam database:

```bash
wget -c ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz
gzip -d ko_list.gz

wget -c ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz
tar zxvf profiles.tar.gz
```

For [seed complementarity](../modules/modules.md), if you do not already have genome-scale metabolic models (GEMs)
of your genomes and reconstruct them with CarveMe, you need a
[Gurobi Web License Service (WLS) license](https://www.gurobi.com/features/academic-wls-license/), or optionally
CPLEX; see [installing a Gurobi license](../installation.md#install-gurobi-license).


### The HPC case

On an HPC, we suggest you first download the images you need. Make sure you have access to Singularity/Apptainer;
if they are available as a [module](https://hpc-wiki.info/hpc/Modules), load it first.

Then choose where to store the images; you may need to ask your admin, or use a folder of your choice, for example
`${HOME}/.singularity/mtg-images`, which is the default of:

```bash
./get_sif.sh
```

To store them somewhere else, use `--sif-dir`:

```bash
./get_sif.sh --sif-dir /opt/sing/mtg-images
```

Once you have the images, set the containerization technology to Singularity or Apptainer, as shown in the
[Nextflow configuration](#nextflow-configuration). You may also need to set environment variables of
Singularity/Apptainer or Nextflow, for example:

    SINGULARITY_TMPDIR
    SINGULARITY_CACHEDIR
    NXF_SINGULARITY_CACHEDIR
    NXF_TEMP

Follow your admin's instructions on how to use Nextflow and Singularity/Apptainer.

Last, run:

```bash
mkdir $HOME/julia_depot $HOME/julia_tmp
```

This is needed for Julia (FlashWeave) to work inside the Singularity/Apptainer containers.


## Nextflow versions

The pipeline follows the strict syntax that Nextflow uses by default from version 25.10;
it was last tested with Nextflow 26.04.


## Workflows

Four workflows, in the `workflows` folder:

1. catalogue pre-calculations (`precalc.nf`)
2. network annotation using a pre-calculation folder (`annotate.nf`)
3. preparation of amplicon datasets for the on-the-fly version (`prep.nf`); see [the preparation steps](../tutorials_otf/prep.md)
4. FlashWeave, FAPROTAX and manta only (`network.nf`), which [the server](./server.md) runs

`main.nf` runs one of them, chosen with `--workflow`. They are built on **modules** (`modules` folder) and
**sub-workflows** (`workflows/subworkflows` folder), which can also run on their own.

The example parameter files read the test data in `tests/data` and `tests/expected`, and write to `tests/results`
(not tracked); `tests/README.md` describes them. Put databases you download, like the KOfam database, in `databases/`
(not tracked either), or anywhere you like and give their paths.

### `precalc`: one-time pre-calculation

Run the `precalc` workflow once, to compute all pairwise relationships between the genomes of your catalogue.
This builds a set of pre-computed results (annotations, complementarities, etc.).

### `annotate`: multiple annotation sessions

Reuse the pre-calculations to annotate multiple networks, from different abundance tables.

```{important}
The species in your abundance table or network need to be among those of the catalogue.
```

### Your own database

The output of `precalc` can also be the database of [your own server](./server.md), which annotates networks sent
from `MGG` with your genomes, by matching their taxa to the closest of your genomes. Export it with the
`microbetag-export-db` command of the package, with the taxonomy of each genome:

```bash
microbetag-export-db --precalc <precalc output folder> --genome-taxonomy genome_taxonomy.tsv --outdir db_files
```

The `export_db` module runs it in the pipeline too:

```bash
nextflow run modules/export_db/main.nf -params-file modules/export_db/export_db.yaml
```

[Your own server](./server.md) describes how to load these files in a database and serve them.


## Configuration

There are **two** configuration levels, with their own files:

- the Nextflow configuration: `nextflow.config`
- the parameters of a run: a params file (YAML), from the template of its workflow

The microbetag steps run in the `hariszaf/microbetag-nf` image, set by the `microbetag_image` parameter; its version
is the version of the `microbetag` package it holds.
The external tools have their own images (see `containers/` and `get_sif.sh`).

### Nextflow configuration

Make sure that Nextflow uses the containerization technology available on your system, by setting its `enabled`
flag to `true` and the others to `false`. For example, to use Docker:

```groovy
docker.enabled    = true
apptainer.enabled = false
singularity {
    enabled    = false
    autoMounts = true
    cacheDir   = "${System.getenv('HOME')}/.singularity/mtg-images/"
    runOptions = [ ... ]   // binds common HPC folders (tmp, scratch, TMPDIR) when they exist
}
```

```{attention}
In the [HPC case](#the-hpc-case), `cacheDir` needs to be the folder the images were downloaded to
(`--sif-dir` of `get_sif.sh`; by default `${HOME}/.singularity/mtg-images/`).
```

### Parameters

At this level, you choose the steps to run and your input data, in a **params file**. The `microbetag-nf`
script of the `pipeline/` folder writes it for a workflow, with every parameter it takes, its default and its
description, and runs it; it calls `main.nf` for you, from any folder (`microbetag-nf help` lists its commands).

```bash
microbetag-nf write-params precalc my_precalc.yaml   # write the params file of the precalc workflow
                                                     # fill it in: your data, the steps you want
microbetag-nf run my_precalc.yaml                    # run it; Nextflow options (-profile, -resume) go after it
```

The params file sets its workflow (`workflow: precalc`), so it is all the pipeline needs. Without the script:

```bash
nextflow run main.nf --workflow precalc --write-params my_precalc.yaml
nextflow run main.nf -params-file my_precalc.yaml
```

A parameter you leave out takes its default, so a params file may set only what you change; a
group of settings (e.g. `flashweave_args`, `modelseed`) may be given in part too. Before a run, the pipeline
reports parameters it does not know (e.g. a typo) and stops if a required one is missing. To describe a parameter,
or all those of a workflow:

```bash
microbetag-nf describe gurobi_sessions     # or: nextflow run main.nf --help gurobi_sessions
microbetag-nf describe annotate            # or: nextflow run main.nf --help --workflow annotate
```

All of them are listed [below](#all-the-parameters). The params files of the tests, in `tests/params/`, are
examples for the test data.

```{attention}
A switch given on the command line is a string for Nextflow: `--phenotrex false` is `"false"`, which counts as
true. Set switches in the params file.
```

In the `precalc` params file, set the output directory for the pre-calculations, the folder with the genomes
or their protein files (`.faa`), and the number of threads. Then set to `true` or `false` the pre-calculation types:

- genome-based phenotypic predictions (`phenotrex`)
- pathway complementarity (`pathway_compl`)
- seed complementarity (`seed_compl`)

Each comes with intermediate files that you may already have, from previous runs or other tools.
For example, `phenotrex` has the `genotype_file`, with the genome id in the first column and the COG ids found in
the genome after it:

```
#feature_type:eggNOG5-tax-2
D300443:bin_000002.fa   COG1040 COG0329 COG4636 COG1373 COG4461 COG2202 COG2084 COG3600 ...
```

If you have it, give its path; otherwise leave the parameter blank.

For `seed_compl` with CarveMe, give the path to your Gurobi WLS license in the `gurobi_lic` parameter.
A WLS license allows a limited number of sessions at the same time (Gurobi tells how many when they are all
in use: "Too many sessions, ... for a baseline of 2"); set it as `gurobi_sessions`, 2 by default. CarveMe runs
one genome per task, at most `gurobi_sessions` at once, and a task that finds all the sessions in use waits and
tries again (`gurobi_retries`, `gurobi_wait`).


## Run

### .. workflows

With the `precalc` params file filled in (the paths to your genomes, the KOfam database and your Gurobi WLS
license), run the workflow with:

```bash
microbetag-nf run my_precalc.yaml
```

`tests/params/precalc.yaml` runs it on the example genomes in `tests/data/genomes`.

```{attention}
On an HPC, make sure you have access to Nextflow and Singularity/Apptainer in the job too;
e.g. in an SBATCH script, you may need to load their modules before the `nextflow` command.
```

Then annotate the network with:

```bash
microbetag-nf write-params annotate my_annotate.yaml
microbetag-nf run my_annotate.yaml
```

`tests/params/annotate.yaml` is filled in for the example pre-calculations (`tests/expected/precalc`) and the example
abundance table in `tests/data/abd_data/`.
Set `faprotax`, `phen_traits`, `path_compl` and `seed_compl` to `true` or `false`, for the annotation types you want.
For the ones set to `true`, give their parameters too; for example, for `faprotax`, `taxonomy_colname` and
`delimiter`. Otherwise, leave them blank.

To prepare an amplicon dataset for the on-the-fly version, see [the preparation steps](../tutorials_otf/prep.md):

```bash
microbetag-nf write-params prep my_prep.yaml
microbetag-nf run my_prep.yaml
```

### .. sub-workflows

The sub-workflows take the parameters of `precalc` (and its defaults). The params files of their tests,
`tests/params/pathway_compl.yaml` and `tests/params/seed_compl.yaml`, are examples:

```bash
nextflow run workflows/subworkflows/pathway_complementarity/main.nf -params-file tests/params/pathway_compl.yaml
```

```bash
nextflow run workflows/subworkflows/seed_complementarity/main.nf -params-file tests/params/seed_compl.yaml
```

### .. modules

Single modules, in the `modules` folder, run the same way. For example:

```bash
nextflow run modules/faprotax/main.nf -params-file modules/faprotax/faprotax.yaml
```


## All the parameters

The templates, `--help` and this list all come from one file of the pipeline,
[`parameters.yaml`](https://github.com/msysbio/microbetag/blob/develop/pipeline/parameters.yaml). A group of
settings is given as a map, e.g. `flashweave_args.n_obs_min` is:

```yaml
flashweave_args:
  n_obs_min: 16
```

```{include} nf_parameters.md
```
