#!/usr/bin/env nextflow

/* Export the output of the precalc workflow as the files of a microbetag database, to load with
server/services/db/schema/custom_init.sql (see server/services/db/README.md).

Run it after precalc, on its output folder:

nextflow run modules/export_db/main.nf -params-file modules/export_db/export_db.yaml
*/

process EXPORT_DB {

    tag "Export the pre-calculations as the files of a microbetag database"

    publishDir "${params.outdir}", mode: 'copy'
    container params.microbetag_image

    input:
        path precalc, stageAs: 'precalc'
        path genome_taxonomy

    output:
        path "db_files", emit: db_files

    script:
        """
        microbetag-export-db --precalc precalc --genome-taxonomy ${genome_taxonomy} --outdir db_files
        """
}


workflow {
    EXPORT_DB(
        channel.fromPath(params.precalc, type: 'dir', checkIfExists: true),
        channel.fromPath(params.genome_taxonomy, checkIfExists: true)
    )
}
