---
title: microbetag-prep
layout: default
parent: Additional tutorials
nav_order: 1
description: "an example case of how to run the microbetag prep step"
---

# `microbetag` preparation steps

<div style="display: flex; gap: 10px;">
   <a href="https://github.com/hariszaf/microbetag/tree/develop/pipeline" class="btn-green" target="_blank"> pipeline </a>
   <a href="../_static/download/prep/seq_ab_tab.tsv" class="btn-purple"> Tutorial file </a>
</div>


```{important}
**INPUT FILES USED IN THIS TUTORIAL**

 - [`seq_ab_tab.tsv`][1]: our abundance table, which in its last column has the **ASV sequence** instead of a taxonomy
 - [`prep.yaml`][2]: the parameters of the preparation workflow
```

## The *preparation* step

The on-the-fly version of `microbetag` annotates your network on our server, from your abundance table.
Two tasks before the annotation can be computationally expensive, especially for large datasets:

- the **taxonomy assignment** of your ASVs/OTUs with GTDB taxonomies, so that `microbetag` can match them to the
  GTDB genomes of our server as well as possible. It uses the IDTAXA algorithm of the
  <a href="https://bioconductor.org/packages/DECIPHER/" target="_blank">DECIPHER</a> package and a training set
  built from the 16S rRNA gene sequences of the GTDB genomes
  (based on <a href="https://zenodo.org/records/6655692" target="_blank">these reference sequences</a>);
- the **network inference** with FlashWeave, needed when you have more than 1,000 ASVs/OTUs and no network.

You can run them on your computer or an HPC with the **prep workflow** of the [pipeline](../tutorials_local/nf.md),
and then upload its outputs to MGG.

```{note}
Up to version 1.0.1, these steps were in a stand-alone Docker image, `hariszaf/microbetag_prep`.
That image is still available on Docker Hub, but it is no longer developed; the prep workflow replaces it.
```


## What you need

- [Nextflow](https://www.nextflow.io/docs/latest/install.html), and Docker or Singularity/Apptainer;
  see [installing the pipeline](../installation.md#the-pipeline)
- the pipeline, the `pipeline/` folder of the `microbetag` repository:

  ```bash
  git clone https://github.com/hariszaf/microbetag.git
  cd microbetag/pipeline
  ```

- for the taxonomy assignment, the **GTDB 16S training set** (`gtdb_r207_ssu_idtaxa.RData`): IDTAXA's training set,
  built with DECIPHER's `LearnTaxa` from the
  <a href="https://zenodo.org/records/6655692" target="_blank">16S rRNA gene sequences of the GTDB r207 genomes</a>.
  r207 is the GTDB release of our server's database: the taxonomies need to be of the same release, since the
  server matches taxa to genomes by their names
- your **abundance table**, as a `.tsv` or `.csv` file:
  the sequence identifier in the first column, then the samples, and in the **last** column either the
  **sequence** (to assign taxonomies) or the **taxonomy** (if you only need the network)
- optionally, a **metadata file** for FlashWeave; see [the input files](../tutorials_core/input.md#abundance-table)
  and the <a href="https://github.com/meringlab/FlashWeave.jl" target="_blank">FlashWeave documentation</a>

```{important}
Make sure there are no comment lines on top of your abundance table: it needs to start with a header.
Avoid numbers as sequence identifiers (e.g., 5434); use alphanumerics instead (e.g., ASV_5434).
```


## Parameters

The parameters are set in a YAML file; [`prep.yaml`][2] is a template.

|**Parameter**             |**Description**                                                                                         |
|--------------------------|--------------------------------------------------------------------------------------------------------|
|`abundance_file`          | The abundance table |
|`outdir`                  | The output folder |
|`assign_taxonomy`         | Assign GTDB taxonomies to the sequences with IDTAXA (`true`/`false`) |
|`gtdb_16s_training_set`   | The GTDB 16S training set (`gtdb_r207_ssu_idtaxa.RData`) |
|`idtaxa_threshold`        | IDTAXA confidence threshold: 60 (cautious) or 50 (sensible) |
|`build_network`           | Infer the co-occurrence network with FlashWeave (`true`/`false`) |
|`metadata_file`           | Metadata file for FlashWeave; leave empty if not available |
|`flashweave_args`         | FlashWeave arguments, e.g. `sensitive`, `heterogeneous`, `n_obs_min`; see the [FAQ](../faq.md#when-to-enable-the-sensitive-and-heterogeneous-arguments) |
|`threads`                 | Number of threads |


## An example case

We will use the output of a 16S rRNA analysis with <a href="https://benjjneb.github.io/dada2/" target="_blank">DADA2</a>,
exported as a `.tsv` file, and show how to get to a *microbetag*-annotated network with this table as the only input.

```{tip}
An optimal co-occurrence network is essential from a biological point of view.
We strongly suggest you build your network with FlashWeave or any inference tool you prefer, to address the
idiosyncrasy of your data the best you can. The prep workflow passes the `flashweave_args` to FlashWeave;
see the [FAQ](../faq.md#when-to-enable-the-sensitive-and-heterogeneous-arguments) for its most essential parameters,
and the <a href="https://github.com/meringlab/flashweave.jl" target="_blank">FlashWeave documentation</a> for the rest.
```

[Download][1] the table. It has 996 ASVs and 16 samples; the on-the-fly version can infer networks of up to
1,000 sequences, but amplicon analyses often have several thousands, and then the preparation step is needed.
In its last column there is the ASV sequence instead of a taxonomy, since we want GTDB taxonomies, from the 16S rRNA gene,
to map our ASVs to their closest GTDB genomes.

In the `pipeline` folder, create a folder for the data, and move the table, the training set and the [`prep.yaml`][2]
template there:

```bash
mkdir my_microbetag_prep
mv ~/Downloads/seq_ab_tab.tsv ~/Downloads/prep.yaml /path/to/gtdb_r207_ssu_idtaxa.RData my_microbetag_prep/
```

In `prep.yaml`:

- `abundance_file` is `my_microbetag_prep/seq_ab_tab.tsv`;
- we will both assign the GTDB taxonomies and build the network, so `assign_taxonomy` and `build_network` are `true`;
- we have no metadata file, so `metadata_file` is empty;
- the table has 16 samples, fewer than the 20 FlashWeave needs by default, so `n_obs_min` is `16`.

Then run, from the `pipeline` folder:

```bash
nextflow run main.nf -params-file my_microbetag_prep/prep.yaml
```

Based on your input, this can take from minutes to hours. You will find in `my_microbetag_prep/prep_output`
(the `outdir` of `prep.yaml`):

- `GTDB_tax_assigned_abundance_table.tsv`: your abundance table with the GTDB taxonomy of each ASV/OTU in place of its sequence;
- `flashweave/fw_net.txt`: the network inferred by FlashWeave, as a 3-column edge list: the first two columns give a
  pairwise association of two ASVs/OTUs, and the third one its weight, positive for co-occurrence and negative for depletion.

Upload both to MGG, as described in [annotating your own network](./from_net.md).
`MGG` sends them to our server, which matches the ASVs to their closest GTDB genomes by their new taxonomies.

```{important}
**REMOVE TAXA NOT PRESENT IN THE NETWORK**

Once you have your network, you can remove from your abundance table the taxa that are not in it. This will save you quite some time!
```


[1]:../_static/download/prep/seq_ab_tab.tsv
[2]:../_static/download/prep/prep.yaml
