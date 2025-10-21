#!/usr/bin/env nextflow

/* Nextflow script to perform KEGG ORTHOLOGY annotation using HMMER and kofam_scan

Usage: 

nextflow run kofam/kofam.nf --kofam_config kofam/kofam.config
*/

include { readParamsFile } from '../helpers.nf'

// Only read default YAML if user didn't specify a params-file
if (!workflow.commandLine.contains('-params-file')) {
    println "[INFO] No params-file provided, loading default YAML..."
    def new_params = readParamsFile(params.paramsFile)
    params.putAll(new_params)
} else {
    println "[INFO] Using user-provided params-file, skipping default YAML."
}

println "Parameters after merge: ${params}"


process hmmsearch{

    publishDir { "results/hmmsearch/${faa.baseName}" }, mode: 'copy'
    container "staphb/hmmer:3.4"

    input:
        tuple path(faa), path(hmmsearch_sc), path(ko_list), path(hmm_profiles)

    output:
        path "hmmout"
    
    script:
    """
    echo "Running HMMER search on ${faa} using kofam: ${ko_list}}"
    mkdir -p hmmout
    bash ${hmmsearch_sc} ${faa} ${ko_list} ${hmm_profiles}
    mv *.hmmout hmmout/
    """
}

workflow {

    // Create a channel from input genomes
    def hmmsearch_sc_ch = Channel.fromPath('kofam/kofam.sh')
    def genomes_ch      = Channel.fromPath("${params.faa_dir}/*.faa")
    def ko_list_ch      = Channel.fromPath("${params.kegg_list}")
    def hmm_prof_ch     = Channel.fromPath("${params.hmm_profiles}") // We do not use *.hmm here to get the folder

    // Combine the single channels with all genomes
    def inputs_ch = genomes_ch
        .combine(hmmsearch_sc_ch)
        .combine(ko_list_ch)
        .combine(hmm_prof_ch)

    // Run HMMER search with kofam_scan
    hmmsearch(inputs_ch)

}


