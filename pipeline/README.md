# microbetag-nf

A [Nextflow](https://www.nextflow.io/) pipeline that takes you from your genomes to a
[microbetag](https://github.com/hariszaf/microbetag)-annotated co-occurrence network.

It runs the external tools, each in its own container (Prodigal, kofam/hmmsearch, phenotrex,
CarveMe/gapseq/ModelSEEDpy, FAPROTAX, FlashWeave, manta), and then the `microbetag` commands, which compute
the pathway and seed complementarities and build the annotated network.

**Documentation:** how to set up and run the pipeline is described on the microbetag documentation site,
[microbetag as a Nextflow pipeline](https://microbetag.readthedocs.io/en/latest/tutorials_local/nf.html).
Its source is `docs/tutorials_local/nf.md` in the microbetag repository; please update it there.

Quick start, with Docker and the example genomes in `data/genomes`:

```bash
nextflow run workflows/precalc.nf -params-file workflows/precalc.yaml -entry MICROBETAG_PRECALC
nextflow run workflows/annotate.nf -params-file workflows/annotate.yaml -entry MICROBETAG_ANNOTATE
```
