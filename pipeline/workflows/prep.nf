#!/usr/bin/env nextflow

/* Prepare an amplicon dataset for the on-the-fly version of microbetag (Cytoscape MGG app):

- assign GTDB taxonomies to the sequences with IDTAXA, so that microbetag can match them to the
  GTDB genomes of microbetagDB, and/or
- infer the co-occurrence network with FlashWeave, for datasets too large for the server.

Upload the outputs to MGG.

Usage:

nextflow run main.nf -params-file workflows/prep.yaml --workflow prep
*/

include { TAXONOMY } from '../modules/taxonomy/'
include { FORMAT_FW; RUN_FW; RUN_FW_METADATA } from '../modules/flashweave/'


workflow MICROBETAG_PREP {

    if (!params.assign_taxonomy && !params.build_network) {
        error "Set assign_taxonomy and/or build_network to true."
    }

    abundance_ch = channel.fromPath(params.abundance_file, checkIfExists: true)

    if (params.assign_taxonomy) {
        log.info "Assigning GTDB taxonomies to the sequences with IDTAXA."
        training_set = channel.fromPath(params.gtdb_16s_training_set, checkIfExists: true)
        TAXONOMY(abundance_ch, training_set)
    }

    if (params.build_network) {
        log.info "Inferring the co-occurrence network with FlashWeave."

        // The last column (the sequences) is dropped when formatting the table for FlashWeave
        formatted_table = FORMAT_FW(abundance_ch)

        if (params.metadata_file) {
            RUN_FW_METADATA(formatted_table, channel.fromPath(params.metadata_file, checkIfExists: true))
        } else {
            RUN_FW(formatted_table)
        }
    }
}
