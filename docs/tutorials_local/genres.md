---
title: Genome-scale metabolic models
layout: default
parent: Additional tutorials
nav_order: 3
description: "the GEMs microbetag uses for seed complementarities"
---


# Genome-scale metabolic models (GEMs)

`microbetag` computes seed complementarities from one genome-scale metabolic model per genome.
It does not build them: the [microbetag-nf](./nf.md) pipeline does, with one of

- <a href="https://carveme.readthedocs.io/en/latest/" target="_blank">CarveMe</a>: fast, works from DNA or protein
  sequences, uses BiGG identifiers; needs a <a href="https://www.gurobi.com/academia/" target="_blank">Gurobi</a>
  (or CPLEX) license;
- <a href="https://gapseq.readthedocs.io" target="_blank">gapseq</a>;
- <a href="https://github.com/ModelSEED/ModelSEEDpy" target="_blank">ModelSEEDpy</a>, uses ModelSEED identifiers;
  needs RAST-annotated genomes, and gap-filling with DNNGIOR.

You can also use models you built yourself, or curated ones, for example to match the medium of your experiment.

## What `microbetag` expects

- A folder with **one SBML (`.xml`) file per genome**, named after the sequence id of the genome in the network
  (`bin_101.xml` for the node `bin_101`). Point `genome_scale_models` in the configuration file to it.
- All models in the **same namespace**: BiGG (CarveMe, gapseq) or ModelSEED (ModelSEEDpy).
  `microbetag` recognises it from the compound identifiers and maps BiGG compounds to ModelSEED for the complements.

```yaml
seed_complementarity:
  value: true
genome_scale_models:
  dir_path: gems
```

The seed sets of the models are written to the `seeds_complementarity` folder of the output directory.
To compute the scores again without the models, point `prev_conf` and `prev_nonseeds` to the
`confidenceDic.json` and `nonSeedSetDic.json` files there.
See [.. using custom genomes](./local.md#seeds_complementarity) for the outputs.
