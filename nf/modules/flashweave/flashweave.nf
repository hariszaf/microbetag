#!/usr/bin/env nextflow

/* Nextflow script to run FlashWeave 

Usage: 

nextflow run flashweave/flashweave.nf  --flashweave_config flashweave/flashweave.config
*/

include { READPARAMSFILE } from '../helpers.nf'
import groovy.json.JsonOutput

// Only read default YAML if user didn't specify a params-file
if (!workflow.commandLine.contains('-params-file')) {
    println "[INFO] No params-file provided, loading default YAML..."
    def new_params = READPARAMSFILE(params.paramsFile)
    params.putAll(new_params)
} else {
    println "[INFO] Using user-provided params-file, skipping default YAML."
}


println "params: ${params}"


// ---------------- Processes ----------------

process formatFlashWeave {
    publishDir 'results', mode: 'copy'
    container "hariszaf/microbetag-nf:0.1.0"

    input:
        path abundance_file
        path fw_formt_script

    output:
        file 'flashweave_abd_table.tsv'
        // file "emm"

    script:
    """
    python3 ${fw_formt_script} ${abundance_file} > flashweave_abd_table.tsv 2>emm
    """
}

process runFW {

    publishDir "${params.outdir}", mode: 'copy'
    container "hariszaf/flashweave:0.19.2"

    input:
        file formatted_table
        file fw_run_sc

    output:
        file 'flashweave.edgelist'

    script:
        """
        julia ${fw_run_sc} \
        --input ${formatted_table} \
        --fw_args '${groovy.json.JsonOutput.toJson(params.flashweave_args)}'
        """
}


process runFWMetadata {

    publishDir "${params.outdir}", mode: 'copy'
    container "hariszaf/flashweave:0.19.2"

    input:
        file formatted_table
        file fw_run_sc
        file metadata_file

    output:
        file 'flashweave.edgelist'

    script:
        """
        julia ${fw_run_sc} \
        --input ${formatted_table} \
        --metadata ${metadata_file} \
        --fw_args '${groovy.json.JsonOutput.toJson(params.flashweave_args)}'
        """   
}

// ---------------- Workflow definition ----------------

workflow {

    abundance_ch = Channel.fromPath(params.abundance_file)
    format_sc_ch = Channel.fromPath('modules/flashweave/format.py')
    fw_sc_ch     = Channel.fromPath('modules/flashweave/run_fw.jl')

    // Step 1: format input table
    formatted_table = formatFlashWeave(abundance_ch, format_sc_ch)

    // Step 2: run FlashWeave on formatted table 
    metadata_val = params.metadata_file ?: ""


    if (metadata_val) {
        // metadata_val is non-empty → run process that uses metadata
        metadata_ch = Channel.fromPath(metadata_val)
        runFWMetadata(formatted_table, fw_sc_ch, metadata_ch)
    } else {
        // metadata_val is empty → run process that does not use metadata
        runFW(formatted_table, fw_sc_ch)
    }
}
