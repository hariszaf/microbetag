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


## Repository layout

```
main.nf                  entry point: runs the workflow chosen with --workflow
nextflow.config          container engine (Docker, Singularity, Apptainer) and the microbetag-nf image version
workflows/               the three workflows and their example parameter files
  precalc.nf/.yaml         pre-calculations for a genome catalogue
  annotate.nf/.yaml        annotate a network with a pre-calculation folder
  prep.nf/.yaml            prepare an amplicon dataset for the on-the-fly version of microbetag
  subworkflows/            pathway and seed complementarity, used by precalc.nf; also runnable on their own
modules/                 one folder per step (a process and an example parameter file), runnable on their own
scripts/                 helper scripts, by the image they run in:
  microbetag-nf/           copied into the microbetag-nf image (FlashWeave input, KO table, taxonomy table)
  flashweave/              copied into the FlashWeave image (runs FlashWeave)
  taxonomy/                IDTAXA, passed to the DECIPHER image
containers/              Dockerfiles of the images built for the pipeline
envs/                    conda environments of the phenotrex, ModelSEEDpy and DNNGIOR images
assets/                  FAPROTAX 1.2.10
get_sif.sh               download all images as Singularity/Apptainer images, e.g. on an HPC
tests/                   example data, reference outputs and test runs (below)
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
| `nextflow run main.nf -params-file workflows/annotate.yaml --workflow annotate` | `tests/data/abd_data/abund.tsv` and the pre-calculations in `tests/expected/precalc` | `tests/results/annotate` | `tests/expected/annotate` |
| `nextflow run main.nf -params-file workflows/precalc.yaml --workflow precalc` | `tests/data/genomes`; needs the KOfam database and, for CarveMe, a Gurobi license | `tests/results/precalc` | `tests/expected/precalc` |
| `nextflow run main.nf -params-file workflows/prep.yaml --workflow prep` | `tests/data/prep`; needs the GTDB 16S training set in `databases/` | `tests/results/prep` | `tests/data/prep/expected_GTDB_tax_assigned_abundance_table.tsv` |
| `nextflow run workflows/subworkflows/<name>/main.nf -params-file workflows/subworkflows/<name>/<pc or sc>.yaml` | the inputs named in its `.yaml` | `tests/results/subworkflows` | the matching folder of `tests/expected/precalc` |
| `nextflow run modules/<step>/main.nf -params-file modules/<step>/<step>.yaml` | the inputs named in its `.yaml` | `tests/results/modules` | the matching folder of `tests/expected/precalc` |

The `annotate` example is the quickest check that the pipeline works (a few minutes): it infers the network with
FlashWeave, runs FAPROTAX, and annotates the network with the example pre-calculations.

```bash
nextflow run main.nf -params-file workflows/annotate.yaml --workflow annotate

# the inferred network should have the same edges as the reference one
# (the weights may differ in their last decimals)
diff <(cut -f1,2 tests/results/annotate/flashweave/fw_net.txt | sort) <(cut -f1,2 tests/expected/annotate/flashweave/fw_net.txt | sort)
```

`tests/README.md` and `tests/data/README.md` describe the test data in more detail.
For your own data, copy a parameter file, point it to your files and set your own output folder.
