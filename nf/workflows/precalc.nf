#!/usr/bin/env nextflow

/* Nextflow script to run microbetag precalculations for a genome catalogue

Usage: 

nextflow run precalc.nf  -params-file params/precalc.yaml

*/

include { readParamsFile } from '../modules/helpers.nf'
include { PHENOTREX } from '../modules/phenotrex/phenotrex.nf'
include { PATHWAY_COMPLEMENTARITY } from '../subworkflows/pathway_complementarity/main'
include { SEED_COMPLEMENTARITY } from '../subworkflows/seed_complementarity/main'

import groovy.json.JsonOutput

// Only read default YAML if user didn't specify a params-file
if (!workflow.commandLine.contains('-params-file')) {
    def new_params = readParamsFile(params.paramsFile)
    params.putAll(new_params)
}


workflow MICROBETAG_PRECALC {


    if ( params.phenotrex ) {

        PHENOTREX()

    }


    if (params.pathway_compl ) {

        PATHWAY_COMPLEMENTARITY()
    }

    if ( params.seed_compl)

    	SEED_COMPLEMENTARITY()
}
