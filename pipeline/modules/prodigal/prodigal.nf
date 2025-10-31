#!/usr/bin/env nextflow

/* Nextflow script to run Prodigal

Usage: 

nextflow run prodigal/prodigal.nf --prodigal_config prodigal/prodigal.config
*/


include { READPARAMSFILE } from '../helpers.nf'

// Only read default YAML if user didn't specify a params-file
if (!workflow.commandLine.contains('-params-file')) {
    println "[INFO] No params-file provided, loading default YAML..."
    def new_params = READPARAMSFILE(params.paramsFile)
    params.putAll(new_params)
} else {
    println "[INFO] Using user-provided params-file, skipping default YAML."
}


process annotate {

    publishDir 'results/annotations', mode: 'copy'
    container "biocontainers/prodigal:v1-2.6.3-4-deb_cv1"
    containerOptions = '-u $(id -u):$(id -g)'

    input:
        path genome

    output:
        path "*.faa"
        path "*.ffn"
        path "*.gbk"


    script:
        def genome_name = genome.name.replaceFirst(/\.[^.]+$/, '')
        """
        prodigal -i ${genome} -p meta -a ${genome_name}.faa -d ${genome_name}.ffn -o ${genome_name}.gbk
        """
}

workflow {

    // Create a channel from input genomes
    def genomes_ch = Channel.fromPath("${params.genomes}/*")

    // Run annotation process
    annotate(genomes_ch)
}
