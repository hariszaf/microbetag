#!/usr/bin/env nextflow

/* Nextflow script to perform Pathway complementarity precalculations 

Usage: 

nextflow run pathway_compl/pathway_compl.nf -c pathway_compl/kofam.config

or

nextflow run pathway_compl/pathway_compl.nf -params-file params/pathway_compl.yaml
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


params.ko_output_file = params.ko_output_file ?: null

if( !params.ko_output_file ) {
    exit 1, "❌ The parameter 'ko_output_file' is required. Example: nextflow run main.nf --ko_output_file myfile.txt"
}



process pathway_compl_precalc {
    tag "Pathway complementarity precalculations"

    publishDir "${params.outdir}/pathway_compl", mode: 'copy'
    container "microbetag"

    input:
    path ko_merged_ch
    path path_compl_sc_ch

    output:
    path params.compl_file
    path params.alts_file

    script:
    """
    python ${path_compl_sc_ch} ${params.alts_file} ${params.compl_file} ${params.tinyurl} ${ko_merged_ch}
    """
}


workflow {

    def ko_merged_ch = Channel.fromPath(params.ko_output_file)
    def path_compl_sc_ch = Channel.fromPath("pathway_compl/pathway_compl.py")

    pathway_compl_precalc(ko_merged_ch, path_compl_sc_ch)


}


