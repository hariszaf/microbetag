#!/usr/bin/env nextflow

/* Nextflow script to microbetag-annotate a network

Usage: 

nextflow run modules/mtg_annotate/main.nf -params-file modules/mtg_annotate/mtg_annotate.yaml
*/


// process ANNOTATE_NETWORK {

//     tag "Build microbetag-annotated network (.cx2 file) "

//     publishDir "${params.outdir}", mode: 'copy'
//     container params.microbetag_image

//     input:
//         path yaml
//         path network
//         path precalc
//         path inDir
//         val fapro_tables
//         val net_cluster

//     output:
//         path "*.cx2", emit: mtg_net

//     script:
//         """
//         cmd="mtg_annotate.py --config_file ${yaml} --network ${network}"

//         [[ -s "${fapro_tables}" ]] && cmd+=" --faprotax ${fapro_tables}"
//         [[ -s "${net_cluster}" ]] && cmd+=" --clustered ${net_cluster}"

//         echo "Running: \$cmd"
//         eval "\$cmd"
//         """
// }
process ANNOTATE_NETWORK {

    tag "Build microbetag-annotated network (.cx2 file)"
    publishDir "${params.outdir}", mode: 'copy'
    container params.microbetag_image

    input:
        path yaml
        path network
        path precalc
        path inDir
        val fapro_tables
        val net_cluster

    output:
        path "*.cx2", emit: mtg_net

    script:
        """

        # Build the command
        cmd="microbetag-annotate --config_file ${yaml} --network ${network}"

        [[ -s "${fapro_tables}" ]] && cmd+=" --faprotax ${fapro_tables}"
        [[ -s "${net_cluster}" ]] && cmd+=" --clustered ${net_cluster}"

        echo "Running: \$cmd"
        eval "\$cmd"
        """
}


workflow {

    // The params file is passed on to microbetag, which reads the annotation settings from it
    def cmdTokens   = workflow.commandLine.tokenize()
    def paramsIndex = cmdTokens.indexOf('-params-file')
    if (paramsIndex == -1) {
        error "Provide the parameters with -params-file, e.g. -params-file modules/mtg_annotate/mtg_annotate.yaml"
    }

    ANNOTATE_NETWORK(
        channel.fromPath(cmdTokens[paramsIndex + 1], checkIfExists: true),
        channel.fromPath(params.network, checkIfExists: true),
        channel.fromPath(params.precalculations, checkIfExists: true),
        channel.fromPath(params.indir, checkIfExists: true),
        params.faprotax_sub_tables ?: "",
        params.manta_net ?: ""
    )
}



// 
        // # --> change paths to be relative to indir

        // # Extract indir value from YAML
        // indir=\$(grep "^indir:" ${yaml} | awk '{print \$2}')

        // # Get last directory name
        // base=\$(basename "\$indir")

        // # Rewrite the YAML file with corrected paths
        // awk -v base="\$base" '
        //     \$1=="indir:" {
        //         print "indir: " base
        //         next
        //     }

        //     \$1=="abundance_file:" {
        //         # extract the path after the key
        //         path=\$2
        //         sub(".*" base "/", base "/", path)
        //         print "abundance_file: " path
        //         next
        //     }

        //     \$1=="metadata_file:" {
        //         if (\$2 != "") {
        //             path=\$2
        //             sub(".*" base "/", base "/", path)
        //             print "metadata_file: " path
        //         } else {
        //             print
        //         }
        //         next
        //     }

        //     \$1=="network:" {
        //         if (\$2 != "") {
        //             path=\$2
        //             sub(".*" base "/", base "/", path)
        //             print "network: " path
        //         } else {
        //             print
        //         }
        //         next
        //     }

        //     { print }
        // ' ${yaml} > tmp.yaml

        // mv tmp.yaml ${yaml}
