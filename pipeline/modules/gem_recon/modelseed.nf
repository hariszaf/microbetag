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


// RAST, on its public server: is it reachable, before the genomes are sent?
process MODELSEED_RAST_CHECK {

    tag "Is RAST reachable?"

    container "hariszaf/modelseed:0.4.2"

    input:
    path faa

    output:
    val true, emit: ok

    script:
    """
    rast_annotate.py --faa ${faa} --check --retries 3 || {
        echo "RAST (https://tutorial.theseed.org/services/genome_annotation) cannot be reached: it is needed to"
        echo "reconstruct GEMs with ModelSEEDpy. Check the internet connection of the compute nodes, or try later."
        exit 1
    } >&2
    """
}


process MODELSEED_RAST {

    tag "RAST annotation of ${base}"

    container "hariszaf/modelseed:0.4.2"
    publishDir path: { "${params.outdir}/reconstructions/${base}" }, mode: 'copy', overwrite: true

    // Not too many genomes at the same time on the public RAST server. A genome RAST failed on (exit status 1),
    // after the retries of rast_annotate.py, is tried again by Nextflow; then, as one with input that is not right
    // (exit status 2: e.g. nucleotides), it stops the run, or is left out with rast.ignore_failed
    maxForks params.modelseed?.rast?.max_forks ?: 4
    errorStrategy {
        task.exitStatus == 1 && task.attempt <= 2 ? 'retry' : (params.modelseed?.rast?.ignore_failed ? 'ignore' : 'terminate')
    }
    maxRetries 2

    input:
    tuple val(base), path(faa)
    val rast_ok

    output:
    tuple val(base), path("${base}-rast.tsv"), emit: rast

    script:
    """
    rast_annotate.py --faa ${faa} --out ${base}-rast.tsv --batch-size ${params.modelseed?.rast?.batch_size ?: 2000}
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
    def rast_todo = faa_ch.filter { base, _faa -> !hasPrev(base, "-rast.tsv") && !hasPrev(base, "-draft.xml") }
    def rast_ok = MODELSEED_RAST_CHECK(rast_todo.first().map { _base, faa -> faa }).ok
    def rast_ch = rast_prev_ch.concat(MODELSEED_RAST(rast_todo, rast_ok).rast)

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
