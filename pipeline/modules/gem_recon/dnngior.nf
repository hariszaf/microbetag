#!/usr/bin/env nextflow

/* Gap-filling of draft GEMs with DNNGIOR (gapfill_with: dnngior), instead of the reconstruction tool's own:
   for the drafts of CarveMe (BiGG), gapseq and ModelSEEDpy (ModelSEED). scripts/dnngior/gapfill.py

DNNGIOR solves with Gurobi: like CarveMe, it needs the WLS licence (params.gurobi_lic), runs at most
params.gurobi_sessions tasks at once, and waits when all the sessions of the licence are in use.
*/

include { gurobiBind; gurobiRetry } from '../helpers.nf'


process DNNGIOR_GAPFILL {

    tag "DNNGIOR gap-filling of ${base}"

    container "hariszaf/dnngior:0.0.3-develop"
    containerOptions { gurobiBind() }
    publishDir "${params.outdir}/reconstructions", mode: 'copy', overwrite: true

    maxForks params.gurobi_sessions ?: 2

    input:
    // The draft is staged in its own folder: the gap-filled model gets the name a draft may have (CarveMe)
    tuple val(base), path(draft, stageAs: 'draft/*')
    // Optional: an empty list for a medium with every compound
    path medium

    output:
    path "${base}.xml", emit: model

    script:
    def medium_arg = medium ? "--medium ${medium}" : ""
    """
    export GRB_LICENSE_FILE=/opt/gurobi/gurobi.lic
    ${gurobiRetry()}
    gurobi_run gapfill.py --model ${draft} ${medium_arg} --out ${base}.xml
    """
}


workflow DNNGIOR {

    take:
    // Draft models (SBML); their genome is the file name, without .xml and -draft
    drafts

    main:
    def med_file = params.dnngior?.medium?.file ?: null
    def med_ch   = med_file ? channel.value(file(med_file, checkIfExists: true)) : channel.value([])
    def models   = DNNGIOR_GAPFILL(
        drafts.map { draft -> [draft.getBaseName().replaceFirst(/-draft$/, ''), draft] },
        med_ch
    )

    emit:
    model = models.model
}
