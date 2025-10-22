#!/usr/bin/env nextflow

/* Nextflow script to run FlashWeave 

Usage: 

nextflow run flashweave/flashweave.nf  --flashweave_config flashweave/flashweave.config
*/

def flashweave_params = [:]
params.flashweave_config = "flashweave/flashweave.config"  // default path

// read the config file
new File(params.flashweave_config).eachLine { line ->
    line = line.replaceAll(/\r/, '').trim()  // remove CRs and spaces
    if (!line.startsWith('#') && line) {
        def (key, value) = line.split('=').collect{ it.trim() }
        flashweave_params[key] = value
    } 
}

println "flashweave_params: ${flashweave_params}"


// ---------------- Processes ----------------

process formatFlashWeave {
    publishDir 'results', mode: 'copy'
    container "microbetag"

    input:
        path abundance_file
        path fw_formt_script

    output:
        file 'flashweave_abd_table.tsv'

    script:
    """
    python3 ${fw_formt_script} ${abundance_file} > flashweave_abd_table.tsv 2>/dev/null
    """
}

process runFW {

    publishDir 'results', mode: 'copy'
    container 'flashweave'

    input:
        file formatted_table  // <- uses the output from previous process
        file fw_run_sc

    output:
        file 'flashweave.edgelist'

    script:
        """
        julia ${fw_run_sc} \
        --input ${formatted_table} \
        --sensitive ${flashweave_params.sensitive} \
        --heterogeneous ${flashweave_params.heterogeneous} \
        --max_k ${flashweave_params.max_k} \
        --n_obs_min ${flashweave_params.n_obs_min} \
        --alpha ${flashweave_params.alpha}
        """
}


process runFWMetadata {

    publishDir 'results', mode: 'copy'
    container 'flashweave'

    input:
        file formatted_table  // <- uses the output from previous process
        file fw_run_sc
        file metadata_file

    output:
        file 'flashweave.edgelist'

    script:
        """
        julia ${fw_run_sc} \
        --input ${formatted_table} \
        --metadata ${metadata_file} \
        --sensitive ${flashweave_params.sensitive} \
        --heterogeneous ${flashweave_params.heterogeneous} \
        --max_k ${flashweave_params.max_k} \
        --n_obs_min ${flashweave_params.n_obs_min} \
        --alpha ${flashweave_params.alpha}
        """   
}

// ---------------- Workflow definition ----------------

workflow {

    abundance_ch = Channel.fromPath(flashweave_params.abundance_file)
    format_sc_ch = Channel.fromPath('flashweave/format.py')
    fw_sc_ch     = Channel.fromPath('flashweave/run_fw.jl')

    // Step 1: format input table
    formatted_table = formatFlashWeave(abundance_ch, format_sc_ch)

    // Step 2: run FlashWeave on formatted table 
    metadata_val = flashweave_params.metadata_file ?: ""

    if (metadata_val) {
        // metadata_val is non-empty → run process that uses metadata
        metadata_ch = Channel.fromPath(metadata_val)
        runFWMetadata(formatted_table, fw_sc_ch, metadata_ch)
    } else {
        // metadata_val is empty → run process that does not use metadata
        runFW(formatted_table, fw_sc_ch)
    }
}
