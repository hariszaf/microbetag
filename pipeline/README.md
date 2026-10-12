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

Quick start, with Docker and the example data (tested with Nextflow 26.04):

```bash
nextflow run main.nf -params-file workflows/precalc.yaml --workflow precalc
nextflow run main.nf -params-file workflows/annotate.yaml --workflow annotate
```

To prepare an amplicon dataset for the on-the-fly version (GTDB taxonomies with IDTAXA, network with FlashWeave):

```bash
nextflow run main.nf -params-file workflows/prep.yaml --workflow prep
```
