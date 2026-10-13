#!/usr/bin/env nextflow

/* Nextflow script to perform seed complementarity precalculations 

Usage: 

nextflow run seed_compl/seed_compl.nf -c seed_compl/seed_compl.config

or

nextflow run seed_compl/seed_compl.nf -params-file params/seed_compl.yaml
*/

include  { getInputFiles; sanitizeChannel; gurobiBind; gurobiRetry; SAFENAME_FILES; GUNZIP } from '../helpers.nf'
include { GAPSEQ } from './gapseq'
include { MODELSEED } from './modelseed'
include { DNNGIOR } from './dnngior'
include { PRODIGAL } from '../prodigal'





process CARVE {

    tag "GEM reconstruction with CarveMe"

    // The models; with gapfill_with: dnngior, the drafts DNNGIOR gap-fills
    publishDir path: { params.gapfill_with == 'dnngior' ? "${params.outdir}/reconstructions/carveme_drafts" : "${params.outdir}/reconstructions" },
        mode: 'copy', overwrite: true
    container "hariszaf/carveme:1.6.6"
    containerOptions { gurobiBind() }

    // One Gurobi session per task: at most as many tasks as the licence has sessions
    maxForks params.gurobi_sessions ?: 2

    input:
        path in_chunk
        val is_faa

    output:
        path "*.xml", emit: carve_model

    when:
        params.recon_with == "carveme"

    script:
        """
        export GRB_LICENSE_FILE=/opt/gurobi/gurobi.lic
        ${gurobiRetry()}

        # Create a clean file list to safely iterate
        for f in ${in_chunk}; do
            echo "\$f" >> chunk.list
        done

        while read f; do
            base=\$(basename "\$f")
            base="\${base%.*}"

            if [ ${is_faa} == "true" ]; then
                echo "Running CarveMe on \$f (FAA mode)" >&2
                gurobi_run carve --solver gurobi --output "\${base}.xml" "\$f" || exit \$?
            else
                echo "Running CarveMe on \$f (DNA mode)" >&2
                gurobi_run carve --dna --solver gurobi --output "\${base}.xml" "\$f" || exit \$?
            fi

        done < chunk.list
        """
}

workflow {

    def recon_in; def prep; def decomp
    def mapped_files; def models

    // sanitize filenames if for example filenames like BATCH:set1.fastq
    recon_in = getInputFiles(params.recon_files)
    prep     = SAFENAME_FILES(recon_in.files_ch)

    // Add a separate mapping step
    mapped_files = prep.map { orig_name, file ->
        def orig_name_decomp = orig_name.endsWith(".gz") ? orig_name[0..-4] : orig_name
        tuple(orig_name, orig_name_decomp, file)
    }

    decomp = GUNZIP(mapped_files)
    decomp_ch = decomp.collect()


    if (params.recon_with == "carveme") {
        // One genome per task; CARVE runs at most gurobi_sessions of them at the same time
        models = CARVE(decomp_ch.flatten(), recon_in.is_faa).carve_model

    } else if (params.recon_with == "gapseq") {

        // gapseq reconstructs from the genomes (DNA), one genome per task (gapseq.nf)
        models = GAPSEQ(decomp_ch.flatten().map { genome -> [genome.getBaseName(), genome] }).model

    } else if (params.recon_with == "modelseed") {

        // ModelSEEDpy reconstructs from proteins: genomes (DNA) are translated with Prodigal first
        def files_ch = decomp_ch.flatten()
        def faa_ch = recon_in.is_faa ? files_ch : PRODIGAL(files_ch).faa

        models = MODELSEED(faa_ch.map { faa -> [faa.getBaseName(), faa] }).model

    } else {
        error "recon_with needs to be carveme, gapseq or modelseed (got: ${params.recon_with})"
    }

    // gapfill_with: dnngior -- the tools gave their drafts, which DNNGIOR gap-fills
    if (params.gapfill_with == 'dnngior') {
        DNNGIOR(models)
    }

}
