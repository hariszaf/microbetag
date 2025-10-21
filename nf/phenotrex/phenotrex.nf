#!/usr/bin/env nextflow

/* Nextflow script to run Phenotrex

Usage: 

nextflow run phenotrex/phenotrex.nf  --phenotrex_config phenotrex/phenotrex.config
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

process genotypes {

    publishDir 'results', mode: 'copy'
    container "phenotrex"

    input:
        path genome_files

    output:
        file 'eggnog.genotype'

    script:
    """
        phenotrex compute-genotype --out eggnog.genotype --threads ${params.threads} ${genome_files.join(' ')}
    """
}


process predict {

    publishDir 'results', mode: 'copy'
    container "phenotrex"

    input:
        tuple path(class_file), path(genotype_file)

    output:
        file "predictions/${class_file.simpleName}.tsv"

    script:
    """
        mkdir -p predictions/
        phenotrex predict --classifier ${class_file} --genotype ${genotype_file} --min_proba ${params.min_proba} --verb > predictions/${class_file.simpleName}.tsv
    """

}


workflow {

    // Step 0: check if user provided precomputed genotype file
    def use_precomputed = params.containsKey('genotype_file') && file(params.genotype_file).exists()
    def classes_ch      = Channel.fromPath("phenotrex/classes/*.pkl")

    // Step 1: generate genotypes from genomes
    if (use_precomputed) {
        println "[INFO] Using user-provided genotype file: ${params.genotype_file}"
        genotype = Channel.fromPath(params.genotype_file)
    } else {
        println "[INFO] No user-provided genotype file found, computing genotypes from genomes."        
        // The .collect() aggregates all files into a single list. The process will receive all files at once.
        genomes_ch = Channel.fromPath("${params.genomes}/*").collect()
        genotype   = genotypes(genomes_ch)

    }

    // Step 2: predict phenotypes from genotypes
    // Each file in the directory becomes one item in the channel. The process will run once per file.
    predict(classes_ch.combine(genotype))
}