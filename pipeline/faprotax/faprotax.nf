#!/usr/bin/env nextflow

/* Nextflow script to run FAPROTAX

Usage: 

nextflow run faprotax/faprotax.nf --faprotax_config faprotax/faprotax.config
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


process faprotax {

    publishDir 'results/faprotax', mode: 'copy'
    container "microbetag"

    input:
        path abundance_table
    
    output:
        path "functional_otu_table.tsv"
        path "sub_tables"
    
    script:
        """
        python /workspace/microbetag/mtg_maps_models/FAPROTAX_1.2.10/collapse_table.py \
            -i ${abundance_table} \
            -o functional_otu_table.tsv \
            -g /workspace/microbetag/mtg_maps_models/FAPROTAX_1.2.10/FAPROTAX.txt \
            --table_delimiter ${params.delimiter} \
            -c "#" \
            -d ${params.taxonomy_colname} \
            -v \
            --force \
            -s sub_tables
        """
}


workflow{

    // Create a channel from input abundance table
    def abundance_ch = Channel.fromPath("${params.abundance_table}")
    faprotax(abundance_ch)

}



