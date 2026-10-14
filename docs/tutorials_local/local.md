---
title: On your own genomes
layout: default
parent: Additional tutorials
nav_order: 3
description: "an example case of how to run microbetag using your own bins/MAGs"
---


# .. using custom genomes

<div style="display: flex; gap: 10px;">
   <a href="https://github.com/hariszaf/microbetag/tree/develop/pipeline" class="btn-green" target="_blank"> pipeline </a>
   <a href="https://github.com/msysbio/microbetag/tree/develop/microbetag/tests/data/microbetag" class="btn-purple" target="_blank"> Tutorial files </a>
</div>


```{note}
This tutorial is for users that have some basic experience working on a terminal.

Contrary to the Cytoscape App tutorials, this scenario is not performed from within Cytoscape:
you run `microbetag` on your computing environment (personal computer, HPC etc.)
and then load the annotated network in Cytoscape.
```

In the [Cytoscape App tutorial](../tutorials_otf/abd_only.md), the taxonomies of the sequences were mapped to
representative GTDB genomes, and `microbetag` used the pre-calculated annotations of those genomes.

With shotgun metagenomics, you may have your own bins or Metagenome-Assembled Genomes (MAGs), or genomes of
the community under study from earlier studies. These can be used directly for the annotation.
This needs much more computing time than our web server can offer, so it runs on your side, in two parts:

1. **[microbetag-nf](./nf.md)**, a Nextflow pipeline, runs the external tools on your genomes:
   gene prediction (Prodigal), KEGG ORTHOLOGY annotation (kofam/hmmsearch), trait prediction (phenotrex),
   genome-scale metabolic model reconstruction (CarveMe, gapseq or ModelSEEDpy), FAPROTAX, and, if needed,
   network inference (FlashWeave) and clustering (manta).
2. **`microbetag`** computes the pathway and seed complementarities from those outputs and builds the
   annotated network.

microbetag-nf runs both parts for you. If you already have some of the tool outputs, for example GEMs or a KEGG
annotation of your genomes from an earlier analysis, you can run `microbetag` on them directly, as below.

If you will annotate many networks with the same genomes, or want others to annotate theirs through `MGG`,
the pre-calculations of the pipeline can also be the database of [your own server](./server.md).


## Input files

The tutorial files are the
<a href="https://github.com/msysbio/microbetag/tree/develop/microbetag/tests/data/microbetag" target="_blank">`tests/data/microbetag`</a>
folder of the `microbetag` repository: 7 bins of the Hessler data set, their network and the tool outputs.

| File | What | Produced by | Used for |
|---|---|---|---|
| `config.yml` | the configuration file; **mandatory** | you | |
| `edgelist.csv` | the co-occurrence network, as a 3-column edge list (node A, node B, weight) | you, or FlashWeave | always |
| `abundance_table.tsv` | sequence ids in the first column, samples, taxonomy in the last column | you | node taxonomy, FAPROTAX |
| `ko_merged.txt` | KEGG ORTHOLOGY terms of each genome | kofam (hmmsearch) | pathway complementarity |
| `gems/` | one genome-scale metabolic model (SBML `.xml`) per genome | CarveMe, gapseq or ModelSEEDpy | seed complementarity |
| `phenotrex_predictions/` | one `.tsv` per trait | phenotrex | phenotypic traits |
| `faprotax_sub_tables/` | one table per process | FAPROTAX | FAPROTAX traits |
| `manta_clustered.cyjs` | the clustered network | manta | network clusters |

Instead of an abundance table, you may give a 2-column file with the sequence id and its taxonomy
(`sequence_id_taxonomy_map`).

```{important}
The genomes are matched to the network nodes by name: the bin ids in `ko_merged.txt`, the names of the
`.xml` files in `gems/` and the identifiers in the phenotrex predictions (without `.fa`) must be the
sequence ids of the network and the abundance table.
```

The formats of the tool outputs:

- **`ko_merged.txt`**: tab-separated, with a header line and one line per annotated gene:

  | bin_id | contig_id | ko_term |
  |:------:|:---------:|:-------:|
  | bin_41 | SCN18_26_2_15_R1_F_scaffold_115_57 | K07586 |
  | bin_48 | SCN18_26_2_15_R4_B_scaffold_93_80  | K08086 |

- **phenotrex predictions** (`phenotrex predict --verb`), e.g. `symbiont.prediction.tsv`:

  ```
  # Trait: Symbiont
  Identifier	Trait present	Confidence
  bin_101.fa	NO	0.643
  bin_45.fa	YES	0.7954
  ```

- **FAPROTAX sub-tables** (`collapse_table.py -s`): one file per process, with the genomes related to it and
  their abundance per sample, e.g. `aerobic_nitrite_oxidation.txt`:

  | record | seqId | sample1 | sample2 | ... |
  |:------:|:-----:|:-------:|:-------:|:---:|
  | d__Bacteria;...;s__Nitrobacter sp001896955 | bin_32 | 43 | 10 | ... |

- **GEMs**: both BiGG (CarveMe, gapseq) and ModelSEED (ModelSEEDpy) models are recognised, from their
  compound identifiers; BiGG compounds are mapped to ModelSEED for the complements.


## The `config.yml` file

`microbetag` takes the paths of these files from a YAML configuration file. Each `microbetag` version comes
with its own template, which `microbetag --write-template config.yml` writes out;
see [the configuration file](./config.md). The one of this tutorial:

```yaml
output_directory:
  dir_path: output

edge_list:
  file_path: inputs/edgelist.csv
abundance_table_file:
  file_path: inputs/abundance_table.tsv

faprotax_annotation:
  value: true
faprotax_sub_tables:
  dir_path: inputs/faprotax_sub_tables

phenotrex_traits:
  value: true
phenotrex_predictions:
  dir_path: inputs/phenotrex_predictions

pathway_complementarity:
  value: true
ko_merged_file:
  file_path: inputs/ko_merged.txt

seed_complementarity:
  value: true
genome_scale_models:
  dir_path: inputs/gems

network_clustering:
  value: true
prev_clustered_network:
  file_path: inputs/manta_clustered.cyjs
```

Paths are relative to the folder of the configuration file.
Set to `false` the annotations you do not have the input for; if an annotation is switched on and its input
is missing, `microbetag` stops and tells you which tool produces it.


## Run

```bash
microbetag --config config.yml
```

On the tutorial files this takes less than a minute on a laptop; the seed complementarities take most of it.
Their time grows with the square of the number of genomes.


## Output

In the `output_directory`:

### The annotated `.cx2` network file

The main output of `microbetag`: the annotated network, `mtag_net_<timestamp>.cx2`, in
[`.cx2`](https://cytoscape.org/cx/cx2/specification/cytoscape-exchange-format-specification-(version-2)/){target="_blank"}
format. [Load it in Cytoscape](../tutorials_core/load.md); after enabling the MGG visual style and the MGG
results panel, you can go through the annotated network.

### `pathway_complementarity`

- `alts.json`: for each genome and KEGG module, the KO terms missing for each alternative path of the module
- `pathCompls.json`: for each pair of genomes, the modules the donor can complete for the beneficiary
- `pathway_complements_extended.json`: the same, with the module descriptions

### `seeds_complementarity`

- `SeedSetDic.json`, `nonSeedSetDic.json`, `confidenceDic.json`: the seed and non-seed sets of each GEM, and
  the confidence of each seed
- `seed_scores.tsv`: the competition and cooperation scores of each pair of genomes
- `seed_complements.pkl`: the seed complements of each pair

To compute the scores again from these sets without the GEMs, point `prev_conf` and `prev_nonseeds` in the
configuration file to `confidenceDic.json` and `nonSeedSetDic.json`.

#### PhyloMInt post process

After computing the seed sets, `microbetag` post-processes the seed and the non-seed sets:
- removing from the seed sets compounds related to environmental metabolites that can be produced in several ways within the cell;
- removing from the non-seed sets compounds that can only be produced by entering the cell from the environment.

The numbers of this post process are recorded in the `log.tsv` file of the `seeds_complementarity` folder:

|model_id | environmental_initial_seeds | non_environmental_initial_seeds | total_initial_seeds | updated_seeds | initial_non_seeds | updated_non_seeds |
|:-------:|:---------------------------:|:-------------------------------:|:-------------------:|:-------------:|:-----------------:|:-----------------:|
| bin_151 |            173              |                57               |            230      |      230      |        879        |       879         |
| bin_38  |            177              |                54               |            231      |      231      |       1211        |      1211         |
| bin_101 |            152              |                51               |            203      |      203      |        845        |       845         |

This post process is needed when the models are gap-filled with a complete medium.

```{note}
One may come with alternative procedures on how to gap fill in terms of minimising the missing potential cross-feedings,
but at the same time not over-predicting such cases.
```

### `parameters.log`

The parameters of the run.

```{important}
`microbetag` reuses the outputs of a previous run with the same output directory.
To compute a step again, for example with different GEMs, delete its output folder first, or use a new output directory.
```
