#!/usr/bin/env nextflow

/* Nextflow script for the network-level steps of microbetag, without the annotation itself:
   network inference with FlashWeave (if no network is given), FAPROTAX and clustering with manta.

   The microbetag server runs it for the on-the-fly version, and then annotates the network with microbetagDB.

Usage:

nextflow run main.nf -params-file workflows/network.yaml --workflow network

Output, in params.outdir: flashweave/fw_net.txt, faprotax/sub_tables/ and manta/manta_net.cyjs
*/

include { FORMAT_FW; RUN_FW; RUN_FW_METADATA } from '../modules/flashweave/'
include { FAPROTAX } from '../modules/faprotax/'
include { PREP_NETWORK; RUN_MANTA } from '../modules/manta/'


workflow MICROBETAG_NETWORK {

    def abundance_ch = channel.fromPath(params.abundance_file, checkIfExists: true)
    def network

    if (params.network) {
        network = channel.fromPath(params.network, checkIfExists: true)
    } else {
        log.info "No network was provided: inferring one with FlashWeave."
        def formatted_table = FORMAT_FW(abundance_ch)
        network = params.metadata_file
            ? RUN_FW_METADATA(formatted_table, channel.fromPath(params.metadata_file, checkIfExists: true))
            : RUN_FW(formatted_table)
    }

    if (params.faprotax) {
        FAPROTAX(abundance_ch)
    }

    if (params.net_cluster) {
        RUN_MANTA(PREP_NETWORK(abundance_ch, network).cyjs_net)
    }
}
