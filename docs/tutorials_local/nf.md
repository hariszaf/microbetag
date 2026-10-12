---
title: microbetag as a Nextflow pipeline
layout: default
parent: Additional tutorials
nav_order: 3
description: "running microbetag on your genomes with Nextflow"
---


# .. as a Nextflow pipeline

<div style="display: flex; gap: 10px;">
   <a href="https://github.com/msysbio/microbetag-nf" class="btn-green" target="_blank"> microbetag-nf </a>
</div>

**microbetag-nf** takes you from your genomes to an annotated network.
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

If you already have the outputs of these tools, you can also run `microbetag` on them directly;
see the tutorial [.. using custom genomes](./local.md).

The pipeline also has a **prep** workflow, which prepares amplicon datasets for the on-the-fly version:
GTDB taxonomies for the sequences, and the network; see [the preparation steps](../tutorials_otf/prep.md).

Paths in this page are relative to the root of the
<a href="https://github.com/msysbio/microbetag-nf" target="_blank">microbetag-nf repository</a>.


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

Three workflows, in the `workflows` folder:

1. catalogue pre-calculations (`precalc.nf`)
2. network annotation using a pre-calculation folder (`annotate.nf`)
3. preparation of amplicon datasets for the on-the-fly version (`prep.nf`); see [the preparation steps](../tutorials_otf/prep.md)

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


## Configuration

There are **two** configuration levels, with their own files:

- the Nextflow configuration: `nextflow.config`
- the `microbetag` configuration: `workflows/precalc.yaml`, `workflows/annotate.yaml` and `workflows/prep.yaml`

The microbetag steps run in the `hariszaf/microbetag-nf` image, set in `params.microbetag_image` of
`nextflow.config`; its version is the version of the `microbetag` package it holds.
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

### `microbetag` configuration

At this level, you choose the steps to run and your input data.

In `workflows/precalc.yaml`, set the output directory for the pre-calculations, the folder with the genomes
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
A WLS license allows up to 5 parallel sessions, so `max_forks` is 5 by default, and the genomes are split
in chunks of `max_forks` for parallel processing.


## Run

### .. workflows

`workflows/precalc.yaml` is filled in for the example genomes in `tests/data/genomes`; set the paths to your genomes,
the KOfam database and your Gurobi WLS license. Run the workflow with:

```bash
nextflow run main.nf -params-file workflows/precalc.yaml --workflow precalc
```

```{attention}
On an HPC, make sure you have access to Nextflow and Singularity/Apptainer in the job too;
e.g. in an SBATCH script, you may need to load their modules before the `nextflow` command.
```

Then annotate the network with:

```bash
nextflow run main.nf -params-file workflows/annotate.yaml --workflow annotate
```

`workflows/annotate.yaml` is filled in for the example pre-calculations (`tests/expected/precalc`) and the example
abundance table in `tests/data/abd_data/`.
Set `faprotax`, `phen_traits`, `path_compl` and `seed_compl` to `true` or `false`, for the annotation types you want.
For the ones set to `true`, give their parameters too; for example, for `faprotax`, `taxonomy_colname` and
`delimiter`. Otherwise, leave them blank.

To prepare an amplicon dataset for the on-the-fly version, see [the preparation steps](../tutorials_otf/prep.md):

```bash
nextflow run main.nf -params-file workflows/prep.yaml --workflow prep
```

### .. sub-workflows

Edit the parameter files of the sub-workflows, `workflows/subworkflows/pathway_complementarity/pc.yaml` and
`workflows/subworkflows/seed_complementarity/sc.yaml`, then run:

```bash
nextflow run workflows/subworkflows/pathway_complementarity/main.nf -params-file workflows/subworkflows/pathway_complementarity/pc.yaml
```

```bash
nextflow run workflows/subworkflows/seed_complementarity/main.nf -params-file workflows/subworkflows/seed_complementarity/sc.yaml
```

### .. modules

Single modules, in the `modules` folder, run the same way. For example:

```bash
nextflow run modules/faprotax/main.nf -params-file modules/faprotax/faprotax.yaml
```
