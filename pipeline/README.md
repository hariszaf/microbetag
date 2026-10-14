# microbetag-nf

A [Nextflow](https://www.nextflow.io/) pipeline that takes you from your genomes to a
[microbetag](https://github.com/hariszaf/microbetag)-annotated co-occurrence network.

It runs the external tools, each in its own container (Prodigal, kofam/hmmsearch, phenotrex,
CarveMe/gapseq/ModelSEEDpy, FAPROTAX, FlashWeave, manta), and then the `microbetag` commands, which compute
the pathway and seed complementarities and build the annotated network. The microbetag steps run in the
`hariszaf/microbetag-nf` image (`containers/Dockerfile.microbetag-nf`), whose version is set in `nextflow.config`.

**Documentation:** how to set up and run the pipeline is described on the microbetag documentation site,
[microbetag as a Nextflow pipeline](https://microbetag.readthedocs.io/en/latest/tutorials_local/nf.html).
Its source is `docs/tutorials_local/nf.md` in the microbetag repository; please update it there.

Tested with Nextflow 26.04 and Docker.


## Requirements

- [Nextflow](https://www.nextflow.io/docs/latest/install.html) and Docker (laptop) or Singularity/Apptainer (HPC)
- the [KOfam database](https://www.genome.jp/ftp/db/kofam/), for the KEGG ORTHOLOGY annotation: `ko_list` and `profiles/`
  (the tests use the small subset in `tests/data/kofam_database`)
- a **Gurobi license**, for the reconstruction of genome-scale metabolic models with CarveMe, and for their
  gap-filling with DNNGIOR (below)


## Gurobi

Two steps of the seed complementarity solve with Gurobi:

- **CarveMe** (`recon_with: carveme`, the default), which reconstructs the models: Gurobi is **required**;
- **DNNGIOR** (`gapfill_with: dnngior`), which gap-fills the drafts of CarveMe, gapseq or ModelSEEDpy.

gapseq (`recon_with: gapseq`) does not use Gurobi: it uses CPLEX if you give its installation (`cplex`), or its
default solver otherwise. ModelSEEDpy does not need Gurobi either.

The pipeline runs CarveMe and DNNGIOR in containers that already have Gurobi (its Python package, `gurobipy`), so
you do **not** need to install Gurobi: you only need a **license file**, given to the pipeline with the
`gurobi_lic` parameter.

**It has to be a Web License Service (WLS) license.** The standard academic license is tied to the machine it was
issued for and does not work inside containers. Academics can get a WLS license for free:

1. Create an account at [gurobi.com](https://www.gurobi.com/) with your institutional email.
2. Request an [academic WLS license](https://www.gurobi.com/features/academic-wls-license/) and download its
   license file from the [Web License Manager](https://license.gurobi.com/); keep it private, it holds your
   license key.
3. Point the pipeline to it, with an absolute path, in your params file (e.g. `precalc.yaml`):

   ```yaml
   gurobi_lic: /home/you/gurobi/gurobi_wls.lic
   ```

   The pipeline mounts it in the CarveMe and DNNGIOR containers as `/opt/gurobi/gurobi.lic`.

A WLS license allows a limited number of sessions at the same time; Gurobi tells how many when they are all in
use ("Too many sessions, ... for a baseline of 2"). Set it as `gurobi_sessions` (2 by default): the tasks that
solve with Gurobi (CarveMe, DNNGIOR) run at most that many at once, one genome per task, and a task that finds
all the sessions in use (e.g. by your other jobs on the same licence) waits and tries again (`gurobi_retries`,
`gurobi_wait`).

## Parameters

A run of the pipeline is described by a **params file** (YAML): the workflow to run and its parameters. You do not
write it from scratch: the pipeline writes a template for the workflow you want, with every parameter it takes,
its default and a description of it, as comments.

1. Write the template of a workflow (here `precalc`) to a file of your choice:

   ```bash
   nextflow run main.nf --workflow precalc --write-params my_precalc.yaml
   ```

2. Open `my_precalc.yaml` and fill it in: the paths to your data (e.g. `genomes`), and the steps you want (e.g.
   `phenotrex: true`). You may delete the parameters you leave at their default, since a parameter missing from
   the file takes its default.

3. Run it. The template already says which workflow it is for (`workflow: precalc`), so the file is all Nextflow
   needs:

   ```bash
   nextflow run main.nf -params-file my_precalc.yaml
   ```

The other workflows work the same way (`--workflow annotate`, `prep` or `network`).

To read the description of a parameter, or of all those of a workflow, without writing a template:
`nextflow run main.nf --help [PARAMETER] [--workflow W]`. All the parameters are in [`parameters.yaml`](./parameters.yaml), where a new parameter is added (with its default,
description and workflows). `python3 tests/check_parameters.py` checks it against the parameters the workflows read
(also in CI).

## Repository layout

```
main.nf                  entry point: runs the workflow of a params file; --write-params, --help
parameters.yaml          every parameter of the workflows: default, description, workflows (the one place)
lib/MtgParameters.groovy reads parameters.yaml: templates, help, checks of a run
nextflow.config          container engine (Docker, Singularity, Apptainer); the defaults, from parameters.yaml
workflows/               the four workflows
  precalc.nf               pre-calculations for a genome catalogue
  annotate.nf              annotate a network with a pre-calculation folder
  prep.nf                  prepare an amplicon dataset for the on-the-fly version of microbetag
  network.nf               FlashWeave, FAPROTAX and manta only; run by the microbetag server (profile `local`)
  subworkflows/            pathway and seed complementarity, used by precalc.nf; also runnable on their own
modules/                 one folder per step (a process and an example parameter file), runnable on their own
scripts/                 helper scripts, by the image they run in:
  microbetag-nf/           copied into the microbetag-nf image (FlashWeave input, KO table, taxonomy table)
  flashweave/              copied into the FlashWeave image (runs FlashWeave)
  taxonomy/                IDTAXA, passed to the DECIPHER image
  modelseed/               copied into the ModelSEEDpy image (RAST annotation, draft reconstruction, gap-filling)
  dnngior/                 copied into the DNNGIOR image (gap-filling of the drafts, with gapfill_with: dnngior)
containers/              Dockerfiles of the images built for the pipeline
envs/                    conda environments of the phenotrex, ModelSEEDpy and DNNGIOR images
assets/                  FAPROTAX 1.2.10
get_sif.sh               download all images as Singularity/Apptainer images, e.g. on an HPC
tests/                   example data, the params files of the tests (params/), reference outputs (below)
databases/               databases you download, e.g. KOfam or the GTDB 16S training set (not tracked)
```


## Tests and examples

The example parameter files are also the tests: each reads from `tests/data` and `tests/expected`, and writes to
`tests/results`, which is not tracked.

```
tests/
  data/        inputs
    genomes/              3 example genomes: the genome catalogue
    abd_data/             an example abundance table (abund.tsv) and metadata for those genomes
    prep/                 an amplicon table with the sequence of each ASV, its metadata, and the output of
                          the former microbetag_prep tool on it
    kofam_database/       KOfam ko_list and a test subset of its profiles (profiles_test.tar.gz; unpack it first)
    manta/                a network and abundance table to test manta
    <tool>/               inputs of single-tool tests (prodigal, phenotrex, carve, modelseed, ...)
  expected/    reference outputs; runs never write here
    precalc/              precalc.nf on tests/data/genomes, one folder per step: prodigal (ORFs),
                          hmmsearch (KO annotation), phenotrex, reconstructions (GEMs), pathway_compl, seed_compl
    annotate/             annotate.nf with expected/precalc and tests/data/abd_data/abund.tsv:
                          flashweave (network), faprotax, and the annotated network (test.cx2)
  results/     where the test runs write (not tracked)
```

| Command | Uses | Writes to | Compare with |
|---|---|---|---|
| `nextflow run main.nf -params-file tests/params/annotate.yaml` | `tests/data/abd_data/abund.tsv` and the pre-calculations in `tests/expected/precalc` | `tests/results/annotate` | `tests/expected/annotate` |
| `nextflow run main.nf -params-file tests/params/precalc.yaml` | `tests/data/genomes`; needs the KOfam database and, for CarveMe, a Gurobi license | `tests/results/precalc` | `tests/expected/precalc` |
| `nextflow run main.nf -params-file tests/params/prep.yaml` | `tests/data/prep`; needs the GTDB 16S training set in `databases/` | `tests/results/prep` | `tests/data/prep/expected_GTDB_tax_assigned_abundance_table.tsv` |
| `nextflow run main.nf -params-file tests/params/network.yaml` | `tests/data/manta/abundance_table.tsv` | `tests/results/network` | (none: FlashWeave's network, FAPROTAX's sub-tables, manta's `.cyjs`) |
| `nextflow run workflows/subworkflows/<name>/main.nf -params-file workflows/subworkflows/<name>/<pc or sc>.yaml` | the inputs named in its `.yaml` | `tests/results/subworkflows` | the matching folder of `tests/expected/precalc` |
| `nextflow run modules/<step>/main.nf -params-file modules/<step>/<step>.yaml` | the inputs named in its `.yaml` | `tests/results/modules` | the matching folder of `tests/expected/precalc` |

The `annotate` example is the quickest check that the pipeline works (a few minutes): it infers the network with
FlashWeave, runs FAPROTAX, and annotates the network with the example pre-calculations.

```bash
nextflow run main.nf -params-file tests/params/annotate.yaml

# the inferred network should have the same edges as the reference one
# (the weights may differ in their last decimals)
diff <(cut -f1,2 tests/results/annotate/flashweave/fw_net.txt | sort) <(cut -f1,2 tests/expected/annotate/flashweave/fw_net.txt | sort)
```

`tests/README.md` and `tests/data/README.md` describe the test data in more detail.
For your own data, copy a parameter file, point it to your files and set your own output folder.
