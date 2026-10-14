#!/usr/bin/env nextflow

/* Prepare an amplicon dataset for the on-the-fly version of microbetag (Cytoscape MGG app):

- assign GTDB taxonomies to the sequences with IDTAXA, so that microbetag can match them to the
  GTDB genomes of microbetagDB, and/or
- infer the co-occurrence network with FlashWeave, for datasets too large for the server.

Upload the outputs to MGG.

Usage:

nextflow run main.nf -params-file tests/params/prep.yaml
*/

include { TAXONOMY } from '../modules/taxonomy/'
include { FORMAT_FW; RUN_FW; RUN_FW_METADATA } from '../modules/flashweave/'


workflow MICROBETAG_PREP {

    if (!params.assign_taxonomy && !params.build_network) {
        error "Set assign_taxonomy and/or build_network to true."
    }

    // Before anything runs: the GTDB 16S reference is not part of the pipeline, and is easy to miss
    if (params.assign_taxonomy && !(params.gtdb_16s_reference && file(params.gtdb_16s_reference).exists())) {
        error """\
            assign_taxonomy needs the GTDB r207 16S reference of IDTAXA (gtdb_16s_reference), not found: ${params.gtdb_16s_reference}
            Download gtdb_r207_ssu_idtaxa.RData from the microbetag data on Zenodo (https://doi.org/10.5281/zenodo.6345382)
            and set its path; or set assign_taxonomy: false to only infer the network.""".stripIndent()
    }

    abundance_ch = channel.fromPath(params.abundance_file, checkIfExists: true)

    if (params.assign_taxonomy) {
        log.info "Assigning GTDB taxonomies to the sequences with IDTAXA."
        reference = channel.fromPath(params.gtdb_16s_reference, checkIfExists: true)
        TAXONOMY(abundance_ch, reference)
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
