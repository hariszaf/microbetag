#!/usr/bin/env nextflow

/* Nextflow script to run FlashWeave 

Usage: 

nextflow run flashweave/flashweave.nf  --flashweave_config flashweave/flashweave.config
*/

include { READPARAMSFILE } from '../helpers.nf'
import groovy.json.JsonOutput

// Only read default YAML if user didn't specify a params-file
if (!workflow.commandLine.contains('-params-file')) {
    def new_params = READPARAMSFILE(params.paramsFile)
    params.putAll(new_params)
}

// ---------------- Processes ----------------

process formatFlashWeave {
    // publishDir '${params.outdir}/flashweave', mode: 'copy'
    container "hariszaf/microbetag-nf:0.1.0"

    input:
        path abundance_file

    output:
        file 'flashweave_abd_table.tsv'

    script:
    """
    format_fw.py ${abundance_file} > flashweave_abd_table.tsv 2>emm
    """
}

process runFW {

    publishDir "${params.outdir}/flashweave", mode: 'copy'
    container "hariszaf/flashweave:0.19.2"

    input:
        file formatted_table

    output:
        file 'flashweave.edgelist'

    script:
        """
        echo "FW ARGS JSON: ${groovy.json.JsonOutput.toJson(params.flashweave_args)}" > emm
        run_fw.jl \
        --input ${formatted_table} \
        --fw_args '${groovy.json.JsonOutput.toJson(params.flashweave_args)}'
        """
}


process runFWMetadata {

    publishDir "${params.outdir}/flashweave", mode: 'copy'
    container "hariszaf/flashweave:0.19.2"

    input:
        file formatted_table
        file metadata_file

    output:
        file 'flashweave.edgelist'

    script:
        """
        println "FW ARGS JSON: ${groovy.json.JsonOutput.toJson(params.flashweave_args)}" > emm

        run_fw.jl \
        --input ${formatted_table} \
        --metadata ${metadata_file} \
        --fw_args '${groovy.json.JsonOutput.toJson(params.flashweave_args)}'
        """   
}

// ---------------- Workflow definition ----------------

workflow {

    abundance_ch = Channel.fromPath(params.abundance_file)

    // Step 1: format input table
    formatted_table = formatFlashWeave(abundance_ch)

    // Step 2: checl if metadata file 
    metadata_val = params.metadata_file ?: ""

    // Step 3: run FlashWeave on formatted table 
    if (metadata_val) {
        // metadata_val is non-empty → run process that uses metadata
        metadata_ch = Channel.fromPath(metadata_val)
        runFWMetadata(formatted_table, metadata_ch)
    } else {
        // metadata_val is empty → run process that does not use metadata
        runFW(formatted_table)
    }
}
