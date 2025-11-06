#!/usr/bin/env nextflow

/* Nextflow script to run Phenotrex

Usage: 

nextflow run phenotrex/main.nf -params-file modules/phenotrex/phenotrex.yaml
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

process GENOTYPE {

    tag "Annotate genomes with EggNOG to build the `genotype` file for phenotrex."

    publishDir "${params.outdir}/phenotrex", mode: 'copy'
    container "hariszaf/phenotrex:0.6.0"

    input:
        path genome_files

    output:
        file 'eggnog.genotype'

    script:
    """
    phenotrex compute-genotype --out eggnog.genotype --threads ${params.threads} ${genome_files.join(' ')}
    """
}


process PREDICT {

    tag "Using the `genotype` file, predict phenotypic traits."

    publishDir "${params.outdir}/phenotrex/predictions", mode: 'copy'
    container "hariszaf/phenotrex:0.6.0"

    input:
        tuple path(class_file), path(genotype_file)

    output:
        file "${class_file.simpleName}.tsv"

    script:
    """
    phenotrex predict --classifier ${class_file} --genotype ${genotype_file} --min_proba ${params.min_proba} --verb > ${class_file.simpleName}.tsv
    """
}


workflow PHENOTREX {

    // Step 0: check if user provided precomputed genotype file
    def use_precomputed = params.containsKey('genotype_file') && params.genotype_file != null && file(params.genotype_file).exists()
    def classes_ch      = Channel.fromPath("modules/phenotrex/classes/*.pkl")

    // Step 1: generate genotypes from genomes
    if (use_precomputed) {
        genotype = Channel.fromPath(params.genotype_file)
    } else {
        // The .collect() aggregates all files into a single list. The process will receive all files at once.
        genomes_ch = Channel.fromPath("${params.genomes}/*").collect()
        genotype   = GENOTYPE(genomes_ch)
    }

    // Step 2: predict phenotypes from genotypes
    // Each file in the directory becomes one item in the channel. The process will run once, per file.
    PREDICT(classes_ch.combine(genotype))
}
