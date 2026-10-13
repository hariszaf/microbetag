---
title: microbetag on containers
layout: default
parent: Additional tutorials
nav_order: 3
description: "running microbetag on a container"
---


# on a container

<div style="display: flex; gap: 10px;">
   <a href="https://hub.docker.com/r/hariszaf/microbetag" class="btn-green" target="_blank"> Docker image </a>
</div>

The `hariszaf/microbetag` image holds the `microbetag` package and its commands. The microbetag-nf pipeline uses
it too; the external tools have images of their own there.

The container needs the **configuration file and the files it points to**, mounted at `/data`.
Paths in the configuration file are relative to its folder, so a configuration that works locally also works
in the container. For example, with the folder of the tutorial [.. using custom genomes](./local.md):

```
my_study/
├── config.yml
└── inputs/
    ├── edgelist.csv
    ├── abundance_table.tsv
    ├── ko_merged.txt
    ├── gems/
    └── ...
```


## .. using Docker

```bash
docker pull hariszaf/microbetag:<version>

docker run --rm -v $PWD/my_study:/data hariszaf/microbetag:<version> microbetag --config /data/config.yml
```

The output directory set in `config.yml` is created inside `my_study`.
Run the image without the `microbetag ...` command to get a shell inside the container.

```{note}
The files the container writes belong to the user the container runs as (`root`, by default).
Add `--user $(id -u):$(id -g)` to `docker run` to have them belong to you.
```


## .. using Singularity/Apptainer

These technologies are widely used on High Performance Computing (HPC) systems.
Convert the Docker image to a Singularity image (`.sif`) once; this does not need `sudo` rights:

```bash
singularity pull microbetag.sif docker://hariszaf/microbetag:<version>
```

Then run `microbetag` with the same mount:

```bash
singularity exec -B $PWD/my_study:/data microbetag.sif microbetag --config /data/config.yml
```
