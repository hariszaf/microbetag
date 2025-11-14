#!/usr/bin/env nextflow

/* Nextflow script to run microbetag precalculations for a genome catalogue

Usage: 

nextflow run precalc.nf  -params-file params/precalc.yaml -entry MICROBETAG_PRECALC

*/

include { readParamsFile } from '../modules/helpers.nf'
include { PHENOTREX } from '../modules/phenotrex/'
include { PRODIGAL } from '../modules/prodigal/'
include { PATHWAY_COMPLEMENTARITY } from '../subworkflows/pathway_complementarity/'
include { SEED_COMPLEMENTARITY } from '../subworkflows/seed_complementarity/'


workflow MICROBETAG_PRECALC {

    if ( params.phenotrex ) {

        PHENOTREX()

    }

    // Both complementarity modules require .faa files at some point. 
    // Make sure either user already provided them from config file, 
    // or microbetag runs Prodigal to get them, before firing the sub-workflows for the complementarities
    def faa_ch
    if (params.pathway_compl || params.seed_compl ) {

        // Handle FAA files
        def has_faa_dir = params.faa_dir && file(params.faa_dir).exists()
        faa_files_empty = true

        if (has_faa_dir) {
            def faa_files_list = file(params.faa_dir).list().findAll { it.endsWith('.faa') }
            faa_files_empty = faa_files_list.isEmpty()
        }

        if (has_faa_dir && !faa_files_empty) {

            // Use existing FAA files
            log.info "✓ Using existing FAA files from: ${params.faa_dir}"
            faa_ch = Channel.fromPath("${params.faa_dir}/*.faa", checkIfExists: true)
        
        } else {

            log.info "○ Running genome annotation to get ORFs with Prodigal "

            def genomes_ch
            genomes_ch = Channel.fromPath("${params.genomes}/*.{fa,fasta}", checkIfExists: true)  // note: do not leave spaces in the regex
            def prodigal_output = PRODIGAL(genomes_ch)
            faa_ch =  prodigal_output.faa

        }
    }

    if (params.pathway_compl ) {

        PATHWAY_COMPLEMENTARITY(faa_ch)
    }

    if ( params.seed_compl) {

        SEED_COMPLEMENTARITY(faa_ch)

    }

}
