#!/usr/bin/env nextflow

/*
microbetag-nf: from your genomes to a microbetag-annotated co-occurrence network.

Choose the workflow with --workflow:

    precalc   pre-calculations for a genome catalogue (annotations, complementarities)
    annotate  annotate a network with a pre-calculation folder
    prep      prepare an amplicon dataset for the on-the-fly version of microbetag
    network   FlashWeave, FAPROTAX and manta only; used by the microbetag server

Usage:

    nextflow run main.nf --workflow precalc --write-params precalc.yaml write the params file (template) of a workflow
    nextflow run main.nf -params-file precalc.yaml                      run it (the file sets its workflow)
    nextflow run main.nf --help [PARAMETER] [--workflow precalc]        describe the parameters

The parameters, their defaults and descriptions are in parameters.yaml (lib/MtgParameters.groovy reads it).
*/

include { MICROBETAG_PRECALC  } from './workflows/precalc'
include { MICROBETAG_ANNOTATE } from './workflows/annotate'
include { MICROBETAG_PREP     } from './workflows/prep'
include { MICROBETAG_NETWORK  } from './workflows/network'


workflow {

    def parameters = new MtgParameters("${projectDir}/parameters.yaml")
    def workflows  = parameters.spec.workflows.keySet()

    if (params.workflow && !(params.workflow in workflows)) {
        error "Choose the workflow to run with --workflow: ${workflows.join(', ')} (got: ${params.workflow})"
    }

    if (params.help) {
        try {
            def name = params.help instanceof String && params.help != 'true' ? params.help : null
            println parameters.describe(name, params.workflow)
        } catch (IllegalArgumentException e) {
            error e.message
        }
        return
    }

    if (params.writeParams) {
        if (!params.workflow) {
            error "Choose the workflow of the params file with --workflow: ${workflows.join(', ')}"
        }
        def target = new File(params.writeParams as String)
        if (target.exists()) {
            error "${target} already exists; choose another file name."
        }
        target.text = parameters.template(params.workflow as String, params.writeParams as String) + '\n'
        println "Params file of the ${params.workflow} workflow written to ${target}; fill it in, then run it with:"
        println "    ${MtgParameters.runCommand(params.writeParams as String)}"
        return
    }

    if (!params.workflow) {
        error "Choose the workflow to run: workflow in the params file (the templates set it), or --workflow: ${workflows.join(', ')}. See: ${MtgParameters.helpCommand()}"
    }

    def unknown = parameters.unknown(params)
    if (unknown) {
        log.warn "Unknown parameters, ignored: ${unknown.join(', ')}. See: ${MtgParameters.helpCommand(params.workflow as String)}"
    }

    parameters.applyDefaults(params)

    def missing = parameters.missing(params, params.workflow as String)
    if (missing) {
        error "The ${params.workflow} workflow needs: ${missing.join(', ')}. See: ${MtgParameters.helpCommand(params.workflow as String)}"
    }

    if (params.workflow == 'precalc') {
        MICROBETAG_PRECALC()
    } else if (params.workflow == 'annotate') {
        MICROBETAG_ANNOTATE()
    } else if (params.workflow == 'prep') {
        MICROBETAG_PREP()
    } else if (params.workflow == 'network') {
        MICROBETAG_NETWORK()
    }
}
