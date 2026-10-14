#!/usr/bin/env nextflow

/* Nextflow script to annotate a co-occurrence network with microbetag, using the pre-calculations of a genome catalogue

Usage:

nextflow run main.nf -params-file annotate.yaml

(nextflow run main.nf --workflow annotate --write-params annotate.yaml writes its params file; tests/params/annotate.yaml is
the example of the tests)

*/

include { FORMAT_FW; RUN_FW; RUN_FW_METADATA } from '../modules/flashweave/'
include { ANNOTATE_NETWORK } from '../modules/mtg_annotate/'
include { FAPROTAX } from '../modules/faprotax/'
include { PREP_NETWORK; RUN_MANTA } from '../modules/manta/'

//MICROBETAG_ANNOTATE
workflow MICROBETAG_ANNOTATE {

    def network
    def net_cluster
    def abundance_ch
    def faprotax_sub_tables

    // microbetag-annotate reads the annotation settings from a YAML file: the parameters of the run, with the
    // defaults of parameters.yaml (main.nf), rather than the params file, which may leave some out
    def params_yaml = file("${workflow.workDir}/annotate_params.yaml")
    params_yaml.text = MtgParameters.toYaml(params) + '\n'
    def yaml = channel.fromPath(params_yaml)

    def abd_data_av = params.abundance_file && file(params.abundance_file).exists()
    def network_av  = params.network && file(params.network).exists()


    if (!network_av && abd_data_av) {

        log.info "No network was provided. microbetag will use FlashWeave to infer one."

        abundance_ch = channel.fromPath(params.abundance_file)  // params.abundance_file

        // Step 1: format input table
        formatted_table = FORMAT_FW(abundance_ch)

        // Step 2: checl if metadata file 
        metadata_val = params.metadata_file ?: ""

        // Step 3: run FlashWeave on formatted table 
        if (metadata_val) {
            // metadata_val is non-empty → run process that uses metadata
            metadata_ch = channel.fromPath(metadata_val)
            network = RUN_FW_METADATA(formatted_table, metadata_ch)
        } else {
            // metadata_val is empty → run process that does not use metadata
            network = RUN_FW(formatted_table)
        }

    } else {

        log.info "A previously inferred network will be used."
        network = channel.fromPath(params.network)
    }


    if (params.faprotax) {

        log.info "Run FAPROTAX against abundance table provided"

        // Create a channel from input abundance table
        abundance_ch = channel.fromPath(params.abundance_file)

        // Run process 
        f = FAPROTAX(abundance_ch)

        faprotax_sub_tables = f.fapro_subtables

    } else {

        faprotax_sub_tables = channel.of([])
    }


    
    if (params.net_cluster) {

        if (params.manta) {

            log.info "A network in .cyjs format was provided and will be used as input for manta."

            f = PREP_NETWORK(abundance_ch, network)

            format_net = f.cyjs_net

            log.info "Run manta..."

            net_cluster = RUN_MANTA(format_net)

        } else {

            log.warn "You have selected clustering your network, but applying your own clustering."
                      "In this case, microbetag expects you to provide your network in .cyjs format."
                      "Otherwise, it will fail when building the annotated network."
        }
    } else {

        net_cluster = channel.of([])
    }


    log.info "Annotating the network..."

    def precalc; precalc = channel.fromPath(params.precalculations)
    def inDir; inDir     = channel.fromPath(params.indir)

    ANNOTATE_NETWORK(yaml, network, precalc, inDir, faprotax_sub_tables, net_cluster)

}

