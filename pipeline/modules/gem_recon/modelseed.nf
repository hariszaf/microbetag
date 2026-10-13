#!/usr/bin/env nextflow

/* GEM reconstruction with ModelSEEDpy, in three steps (the scripts of scripts/modelseed):

    MODELSEED_RAST     RAST annotation of the proteins (needs an internet connection: RAST runs on its server)
    MODELSEED_DRAFT    draft model from the proteins and their RAST annotation, with a ModelSEED template
    MODELSEED_GAPFILL  gap-filling of the draft on a medium (params.modelseed.medium.file), or the complete medium

The steps whose files are in params.modelseed.prev_files (<genome>-rast.tsv, <genome>-draft.xml, in the folder or in
a folder per genome, e.g. the reconstructions folder of a previous run) are skipped.
ModelSEEDpy reconstructs from proteins: genomes (DNA) are translated with Prodigal first.
*/

include { PRODIGAL } from '../prodigal'


// A file of a previous run (params.modelseed.prev_files), e.g. prevFile('bin_101', '-draft.xml'): in the folder
// itself, or in a folder per genome, as this module publishes them (reconstructions/bin_101/bin_101-draft.xml)
def prevFile(base, suffix) {
    def prev = params.modelseed?.prev_files
    if (!prev) {
        return null
    }
    def flat = file("${prev}/${base}${suffix}")
    return flat.exists() ? flat : file("${prev}/${base}/${base}${suffix}")
}

def hasPrev(base, suffix) {
    def f = prevFile(base, suffix)
    return f != null && f.exists()
}


process MODELSEED_RAST {

    tag "RAST annotation of ${base}"

    container "hariszaf/modelseed:0.4.2"
    publishDir path: { "${params.outdir}/reconstructions/${base}" }, mode: 'copy', overwrite: true

    input:
    tuple val(base), path(faa)

    output:
    tuple val(base), path("${base}-rast.tsv"), emit: rast

    script:
    """
    rast_annotate.py --faa ${faa} --out ${base}-rast.tsv
    """
}


process MODELSEED_DRAFT {

    tag "ModelSEEDpy draft reconstruction of ${base}"

    container "hariszaf/modelseed:0.4.2"
    publishDir path: { "${params.outdir}/reconstructions/${base}" }, mode: 'copy', overwrite: true

    input:
    tuple val(base), path(faa), path(rast)

    output:
    tuple val(base), path("${base}-draft.xml"), emit: draft

    script:
    """
    recon.py \
        --faa ${faa} \
        --rast ${rast} \
        --template ${params.modelseed.template} \
        --id ${base} \
        --out ${base}-draft.xml
    """
}


process MODELSEED_GAPFILL {

    tag "ModelSEEDpy gap-filling of ${base}"

    container "hariszaf/modelseed:0.4.2"
    publishDir path: { "${params.outdir}/reconstructions" }, mode: 'copy', overwrite: true

    input:
    tuple val(base), path(draft)
    // Optional: an empty list for the complete medium
    path medium

    output:
    path "${base}.xml", emit: model

    script:
    def medium_arg = medium ? "--medium ${medium}" : ""
    """
    gapfill.py \
        --model ${draft} \
        --template ${params.modelseed.template} \
        ${medium_arg} \
        --out ${base}.xml
    """
}


workflow MODELSEED {

    take:
    // [base, proteins (.faa)] of each genome
    faa_ch

    main:
    // RAST annotation, unless available from a previous run or not needed (draft available)
    def rast_prev_ch = faa_ch
        .filter { base, _faa -> hasPrev(base, "-rast.tsv") }
        .map { base, _faa -> [base, prevFile(base, "-rast.tsv")] }
    def rast_ch = rast_prev_ch.concat(
        MODELSEED_RAST(faa_ch.filter { base, _faa -> !hasPrev(base, "-rast.tsv") && !hasPrev(base, "-draft.xml") }).rast
    )

    // Draft models, unless available from a previous run
    def draft_prev_ch = faa_ch
        .filter { base, _faa -> hasPrev(base, "-draft.xml") }
        .map { base, _faa -> [base, prevFile(base, "-draft.xml")] }
    def draft_input_ch = faa_ch
        .filter { base, _faa -> !hasPrev(base, "-draft.xml") }
        .join(rast_ch)
    def draft_ch = draft_prev_ch.concat(MODELSEED_DRAFT(draft_input_ch).draft)

    // Gap-filling, on the medium if one is given, otherwise on the complete medium
    def med_file = params.modelseed?.medium?.file ?: null
    def med_ch   = med_file ? channel.value(file(med_file, checkIfExists: true)) : channel.value([])
    def models   = MODELSEED_GAPFILL(draft_ch, med_ch)

    emit:
    model = models.model
}


workflow {

    // Proteins (.faa) as they are; genomes (DNA) translated with Prodigal first
    def is_faa = file(params.recon_files).listFiles().any { f -> f.name =~ /\.faa(\.gz)?$/ }
    def faa_files = is_faa
        ? channel.fromPath("${params.recon_files}/*.faa")
        : PRODIGAL(channel.fromPath("${params.recon_files}/*")).faa

    MODELSEED(faa_files.map { faa -> [faa.getBaseName(), faa] })
}
