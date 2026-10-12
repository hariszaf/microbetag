#!/usr/bin/env nextflow

/*
microbetag-nf: from your genomes to a microbetag-annotated co-occurrence network.

Choose the workflow with --workflow:

    precalc   pre-calculations for a genome catalogue (annotations, complementarities)
    annotate  annotate a network with a pre-calculation folder
    prep      prepare an amplicon dataset for the on-the-fly version of microbetag
    network   FlashWeave, FAPROTAX and manta only; used by the microbetag server

Usage:

    nextflow run main.nf -params-file workflows/<workflow>.yaml --workflow <workflow>
*/

include { MICROBETAG_PRECALC  } from './workflows/precalc'
include { MICROBETAG_ANNOTATE } from './workflows/annotate'
include { MICROBETAG_PREP     } from './workflows/prep'
include { MICROBETAG_NETWORK  } from './workflows/network'


workflow {

    if (params.workflow == 'precalc') {
        MICROBETAG_PRECALC()
    } else if (params.workflow == 'annotate') {
        MICROBETAG_ANNOTATE()
    } else if (params.workflow == 'prep') {
        MICROBETAG_PREP()
    } else if (params.workflow == 'network') {
        MICROBETAG_NETWORK()
    } else {
        error "Choose the workflow to run with --workflow: precalc, annotate, prep or network (got: ${params.workflow})"
    }
}
