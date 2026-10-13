---
title: The YAML configuration file
layout: default
parent: Tutorial for running microbetag locally
nav_order: 3
description: >
    microbetag uses a YAML configuration file to provide user's settings 
    In this tutorial we explain when needed and how to fill it in
---


# Configuration file

Every `microbetag` release comes with its YAML configuration template. Get the one of the version you have
installed with:

```bash
microbetag --template config.yml
```

and fill it in.

From version 2.0, the configuration file points to:

- the **network**, as a 3-column edge list (`edge_list`), and an **abundance table** (`abundance_table_file`)
  or a 2-column sequence id to taxonomy file (`sequence_id_taxonomy_map`);
- for each annotation, a switch (`value: true` or `false`) and the **output of the tool** it is computed from:

| Annotation | Switch | Input |
|---|---|---|
| Pathway complementarity | `pathway_complementarity` | `ko_merged_file` |
| Seed complementarity | `seed_complementarity` | `genome_scale_models` |
| Phenotypic traits | `phenotrex_traits` | `phenotrex_predictions` |
| FAPROTAX traits | `faprotax_annotation` | `faprotax_sub_tables` |
| Network clusters | `network_clustering` | `prev_clustered_network` |

Every value is written under a `value`, `file_path` or `dir_path` key, for example:

```yaml
seed_complementarity:
  value: true
genome_scale_models:
  dir_path: gems
```

Paths are relative to the folder of the configuration file, so a configuration file and the files it points to
can be moved together, or mounted together in a container.

The template describes every parameter, including those that let you reuse previous results
(`pc_file`, `prev_conf`, `prev_nonseeds`). The tutorial [.. using custom genomes](./local.md) goes through a
complete example.

## Upgrading from 1.x

Up to version 1.0.4, `microbetag` ran the external tools itself, and the configuration file had their settings
(`bins_fasta`, `kofam_database`, `sc_input_type`, `flashweave_args`, ...). These steps now run in the
[microbetag-nf](./nf.md) pipeline, and each step that used to run a tool takes that tool's output instead:

| Step | 1.0.4 | 2.0 |
|---|---|---|
| Network | inferred with FlashWeave if `edge_list` was missing | `edge_list` is required |
| Pathway complementarity | kofam annotation of `bins_fasta` with `kofam_database` | `ko_merged_file` |
| Seed complementarity | GEMs built from `bins_fasta` (`sc_input_type`, `genre_reconstruction_with`, ...) | `genome_scale_models` |
| Phenotypic traits | phenotrex on `bins_fasta` | `phenotrex_predictions` |
| FAPROTAX | run on the abundance table | `faprotax_sub_tables` |
| Network clustering | manta | `prev_clustered_network` |

The 1.0.4 template also had `pc_percentage`, `prev_calc_path_compl` and `precalulations_only`, which microbetag
did not read; the keys it reads are `pc_percent`, `pc_file` and `precalculations_only`.
