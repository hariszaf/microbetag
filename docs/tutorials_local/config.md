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

A version-specific YAML configuration template accompanies every `microbetag` release, in the
<a href="https://github.com/msysbio/microbetag/tree/develop/config_files" target="_blank">`config_files`</a>
folder of the repository. Use the one of the version you have installed (`microbetag --version`).

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

```{note}
Up to version 1.0.4, `microbetag` ran the external tools itself, and the configuration file had their settings
(`bins_fasta`, `kofam_database`, `sc_input_type`, `flashweave_args`, ...). These steps now run in the
[microbetag-nf](./nf.md) pipeline. The
<a href="https://github.com/msysbio/microbetag/blob/develop/config_files/README.md" target="_blank">`config_files` README</a>
maps the 1.0.4 parameters to their 2.0 inputs.
```
