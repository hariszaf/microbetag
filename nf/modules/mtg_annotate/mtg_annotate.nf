#!/usr/bin/env nextflow

/* Nextflow script to microbetag-annotate a network

Usage: 

nextflow run mtg_annotate/mtg_annotate.nf -c mtg_annotate/mtg_annotate.config
*/

include { readParamsFile } from '../helpers.nf'
import groovy.json.JsonOutput

// Only read default YAML if user didn't specify a params-file
if (!workflow.commandLine.contains('-params-file')) {
    println "[INFO] No params-file provided, loading default YAML..."
    def new_params = readParamsFile(params.paramsFile)
    params.putAll(new_params)
} else {
    println "[INFO] Using user-provided params-file, skipping default YAML."
}


process ANNOTATE_NETWORK {

    tag "Build microbetag-annotated network (.cx2 file) "

    publishDir "${params.outdir}", mode: 'copy'
    container "hariszaf/microbetag-nf:0.1.0"

    input:
    path yaml
    path network
    path precalc
    path inDir
    path fapro_tables

    output:
    path "*.cx2", emit: mtg_net

    script:
    """
        if [ ! -s ${fapro_tables} ]; then
            mtg_annotate.py --config_file ${yaml} --network ${network}
        else
            mtg_annotate.py --config_file ${yaml} --network ${network} --faprotax ${fapro_tables}
        fi
    """
}



workflow {

    def yaml_ch   = Channel.fromPath(params.paramsFile)
    def input_ch  = Channel.fromPath(params.indir)
    def outdir_ch = Channel.fromPath(params.outdir)

    // ANNOTATE_NETWORK(mtg_net_sc_ch, yaml_ch, input_ch, outdir_ch)
    ANNOTATE_NETWORK(yaml_ch, input_ch, outdir_ch)

}
