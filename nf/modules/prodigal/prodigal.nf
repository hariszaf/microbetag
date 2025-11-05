#!/usr/bin/env nextflow

/* Nextflow script to run Prodigal

Usage: 

nextflow run prodigal/prodigal.nf --prodigal_config prodigal/prodigal.config
*/


include { readParamsFile } from '../helpers.nf'

// Only read default YAML if user didn't specify a params-file
if (!workflow.commandLine.contains('-params-file')) {
    def new_params = readParamsFile(params.paramsFile)
    params.putAll(new_params)
}


process PRODIGAL {

    publishDir "${params.outdir}/prodigal", mode: 'copy'
    container "biocontainers/prodigal:v1-2.6.3-4-deb_cv1"
    containerOptions = '-u $(id -u):$(id -g)'

    input:
        path genome

    output:
        path "*.faa", emit: faa
        path "*.ffn", emit: ffn
        path "*.gbk", emit: gbk

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
    // PRODIGAL(genomes_ch)
    faa_ch = PRODIGAL(genomes_ch)

    // debug output in Nextflow log
    faa_ch.view()  

    // Return only the *.faa channel
    return faa_ch

}
