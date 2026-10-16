---
title: Your own microbetag server
layout: default
parent: Additional tutorials
nav_order: 4
description: "running a microbetag server, with a database of your own genomes"
---


# Your own `microbetag` server

<div style="display: flex; gap: 10px;">
   <a href="https://github.com/hariszaf/microbetag/tree/develop/server" class="btn-green" target="_blank"> server </a>
</div>

The **server** is what `MGG` talks to in the [on-the-fly version](../modes.md#mgg-and-a-server). It receives an
abundance table, and optionally a network, and returns the annotated network:

1. it matches each taxon to the closest genomes of its **database**;
2. it runs the `network` workflow of the [pipeline](./nf.md), for what needs to be computed on the data itself:
   FlashWeave, if no network is given, FAPROTAX and manta;
3. it runs the `microbetag` package, which annotates the network with the traits and complementarities
   of the matched genomes, pre-computed in the database.

Our public server has a database of the GTDB representative genomes. You can run a server of your own, for
example for a **genome catalogue of your own**, so that you and your colleagues annotate networks with it, without
running the pipeline for each network. The server is the `server/` folder of the repository.


## What you need

- <a href="https://docs.docker.com/get-docker/" target="_blank">Docker</a>: the server and its database run in containers;
- the repository:

  ```bash
  git clone https://github.com/hariszaf/microbetag.git
  cd microbetag
  ```

- for a database of your own genomes, the [pipeline](./nf.md) and its requirements, to run `precalc` on them.


## 1. Pre-calculations for your genomes

Run the `precalc` workflow of the pipeline on your genomes, as described in [the pipeline](./nf.md#precalc-one-time-pre-calculation):

```bash
cd pipeline
nextflow run main.nf -params-file tests/params/precalc.yaml
```

This is the step that takes time: gene prediction, KEGG annotation, trait prediction and GEM reconstruction of each
genome, and their pairwise complementarities. You run it once.


## 2. The database

Export the output of `precalc` as the files of a database, with the **taxonomy of each genome**: a tab-separated file
with the genome identifier and its lineage, such as the `gtdbtk.*.summary.tsv` of GTDB-Tk, or

```
bin_101	d__Bacteria;p__Firmicutes;c__Clostridia;o__Peptostreptococcales;f__Peptostreptococcaceae;g__Clostridioides;s__Clostridioides difficile
bin_151	d__Bacteria;p__Firmicutes;c__Clostridia;o__Eubacteriales;f__Eubacteriaceae;g__Eubacterium;s__
```

The genome identifiers are the names of the genome files, without their extension (`bin_101.fa` → `bin_101`).
Export with the `microbetag-export-db` command of the package, here from its image:

```bash
docker run --rm -v $PWD:/data -w /data hariszaf/microbetag:<version> \
    microbetag-export-db --precalc <precalc output folder> --genome-taxonomy genome_taxonomy.tsv --outdir db_files
```

`db_files` has one file per table of the database, seed sets included: the server reads the rows of a network's
genomes, and computes their pathway and seed complementarities.

Then start the database: MySQL creates the tables and loads the files at its first start, with
`server/services/db/schema/microbetagDB.sql`:

```bash
docker network create microbetag
docker run -d --name microbetag-db --network microbetag \
    -e MYSQL_ROOT_PASSWORD=<password> \
    -v $PWD/server/services/db/schema/microbetagDB.sql:/docker-entrypoint-initdb.d/microbetagDB.sql:ro \
    -v $PWD/db_files:/var/lib/mysql-files \
    -p 3306:3306 mysql:8.4
```

```{attention}
The `db_files` mount needs to be writable, since MySQL changes its owner when it starts.
On Fedora and other systems with SELinux, add `:Z` to both mounts.
```


## 3. The server

Build its image. It has the `microbetag` package, the tools of the pipeline (FlashWeave, FAPROTAX, manta) and Nextflow,
which runs the pipeline's `network` workflow inside the server's container, without containers of its own.
From the root of the repository:

```bash
docker build -t hariszaf/microbetag:<version> microbetag
docker build -f pipeline/containers/Dockerfile.microbetag-nf -t hariszaf/microbetag-nf:<version> pipeline
docker build -f pipeline/containers/Dockerfile.flashweave -t hariszaf/flashweave:0.19.2 pipeline
docker build -f server/services/web/Dockerfile.v2 -t microbetag-server:<version> .
```

Set its database credentials in a `.env` file, from the template `server/.env.example`:

| Variable | Value |
|---|---|
| `USER_NAME`, `PASSWORD` | the database user and password |
| `HOST`, `PORT` | where the database runs |
| `DB_NAME` | the database, `microbetagDB` with `microbetagDB.sql` |
| `FLASK_SECRET_KEY` | a random string |

and run it. It reaches the database on their shared Docker network
(`microbetag`), so `HOST` is `microbetag-db` and `PORT` is `3306`.

```bash
docker run -d --name microbetag-server --network microbetag -p 5000:5000 --env-file .env \
    -e CATALOGUE=custom \
    microbetag-server:<version>
```

| Variable | Value |
|---|---|
| `CATALOGUE` | `custom`: each taxon is matched to the genomes of the database that share its deepest rank, from its taxonomy; a taxon named after a genome of the database is matched to that genome |
| `CATALOGUE_MIN_RANK` | the shallowest rank a match may share; `genus` by default, so that a taxon that shares only its family with your genomes is not annotated |

To keep the server running after you log out, and start it again at boot, run it with systemd;
`server/README.md` describes how, and how to run it with Docker Compose together with its database.


## 4. Annotate a network

The server annotates the networks it receives at `/upload-abundance-table-dev`, in the format `MGG` sends them:
the abundance table, with the taxonomy in its last column, an optional network, and the parameters.

```bash
curl -X POST -H "Content-Type: application/json" --data @request.json \
    http://localhost:5000/upload-abundance-table-dev -o annotated.cx2
```

```json
{
  "input_type": "network",
  "abundance_data": [["seqId", "S1", "S2", "S3", "taxonomy"],
                     ["bin_101", "12", "0", "31", "Bacteria;Firmicutes;Clostridia;..."],
                     ["ASV_7", "4", "18", "2", "Bacteria;Firmicutes;Clostridia;...;Eubacterium;Eubacterium limosum"]],
  "network": [["node_A", "node_B", "weight"], ["bin_101", "ASV_7", "0.4"]],
  "metadata": null,
  "parameters": ["input_type:network", "taxonomy:other", "delimiter:;",
                 "phen_traits:true", "pathway_complementarity:true", "seed_complementarity:true",
                 "faprotax:true", "network_clustering:false"]
}
```

Without a `network` (`"input_type": "abundance table"`), the server infers it with FlashWeave.
[Load](../tutorials_core/load.md) `annotated.cx2` in Cytoscape and go through it with `MGG`.
