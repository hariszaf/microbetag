#!/usr/bin/env nextflow

/* Nextflow script to run microbetag precalculations for a genome catalogue

Usage: 

nextflow run precalc.nf  -params-file params/precalc.yaml

*/

include { READPARAMSFILE } from './modules/helpers.nf'
include { PHENOTREX } from './modules/phenotrex/phenotrex.nf'
include { KO_ANNOTATE } from './modules/kofam/kofam.nf'
include { ORFS } from './modules/prodigal/prodigal.nf'

import groovy.json.JsonOutput

// Only read default YAML if user didn't specify a params-file
if (!workflow.commandLine.contains('-params-file')) {
    def new_params = READPARAMSFILE(params.paramsFile)
    params.putAll(new_params)
}






workflow {

    if (params.phenotrex) {
        PHENOTREX()
    }


    if (params.pathway_compl && !(params.ko_merged?.trim() && file(params.ko_merged).exists())) {

        KO_ANNOTATE()
    }





}
