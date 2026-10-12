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
        // Optional: an empty list when not available
        path fapro_tables, stageAs: 'faprotax_sub_tables'
        path net_cluster

    output:
        path "*.cx2", emit: mtg_net

    script:
        """
        # The paths in the params file are relative to the launch folder (e.g. tests/expected/precalc/...),
        # while Nextflow stages the input folders under their last name only (precalc/, data/).
        # Link them under their relative paths too, so that microbetag finds the files the params file names.
        link_as() {
            local target="\${1%/}" staged="\$2"
            case "\$target" in /*|"") return ;; esac
            if [ "\$target" != "\$staged" ] && [ ! -e "\$target" ]; then
                mkdir -p "\$(dirname "\$target")"
                ln -s "\$PWD/\$staged" "\$target"
            fi
        }
        link_as "${params.precalculations}" "${precalc}"
        link_as "${params.indir}" "${inDir}"

        # Build the command
        cmd="microbetag-annotate --config_file ${yaml} --network ${network}"

        [[ -n "${fapro_tables}" ]] && cmd+=" --faprotax ${fapro_tables}"
        [[ -n "${net_cluster}" ]] && cmd+=" --clustered ${net_cluster}"

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
        params.faprotax_sub_tables ? channel.fromPath(params.faprotax_sub_tables, checkIfExists: true) : channel.of([]),
        params.manta_net ? channel.fromPath(params.manta_net, checkIfExists: true) : channel.of([])
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
