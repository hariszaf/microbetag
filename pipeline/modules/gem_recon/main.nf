#!/usr/bin/env nextflow

/* Nextflow script to perform seed complementarity precalculations 

Usage: 

nextflow run seed_compl/seed_compl.nf -c seed_compl/seed_compl.config

or

nextflow run seed_compl/seed_compl.nf -params-file params/seed_compl.yaml
*/

include  { getInputFiles; sanitizeChannel; chunkFiles; SAFENAME_FILES; GUNZIP } from '../helpers.nf'
include { GAPSEQ } from './gapseq'
include { MODELSEED } from './modelseed'
include { PRODIGAL } from '../prodigal'


// Mount the Gurobi license in the CarveMe containers
def gurobiBind() {
    return workflow.containerEngine == 'singularity'
        ? "-B ${params.gurobi_lic}:/opt/gurobi/gurobi.lic"
        : "-v ${params.gurobi_lic}:/opt/gurobi/gurobi.lic"
}



process CARVE {

    tag "GEM reconstruction with CarveMe"

    publishDir "${params.outdir}/reconstructions", mode: 'copy', overwrite: true
    container "hariszaf/carveme:1.6.6"
    containerOptions { gurobiBind() }

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

        # Create a clean file list to safely iterate
        for f in ${in_chunk}; do
            echo "\$f" >> chunk.list
        done

        while read f; do
            base=\$(basename "\$f")
            base="\${base%.*}"

            if [ ${is_faa} == "true" ]; then
                echo "Running CarveMe on \$f (FAA mode)" >&2
                carve --solver gurobi --output "\${base}.xml" "\$f"
            else
                echo "Running CarveMe on \$f (DNA mode)" >&2
                carve --dna --solver gurobi --output "\${base}.xml" "\$f"
            fi

        done < chunk.list
        """
}

workflow {

    def recon_in; def prep; def decomp
    def mapped_files;  def infiles_chunk_ch

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
        // Use the COLLECTED channel for chunking, not the original decomp channel
        infiles_chunk_ch = chunkFiles(decomp_ch, params.max_forks)

        // Pass chunks to your carve process
        CARVE(infiles_chunk_ch, recon_in.is_faa)

    } else if (params.recon_with == "gapseq") {

        // gapseq reconstructs from the genomes (DNA), one genome per task (gapseq.nf)
        GAPSEQ(decomp_ch.flatten().map { genome -> [genome.getBaseName(), genome] })

    } else if (params.recon_with == "modelseed") {

        // ModelSEEDpy reconstructs from proteins: genomes (DNA) are translated with Prodigal first
        def files_ch = decomp_ch.flatten()
        def faa_ch = recon_in.is_faa ? files_ch : PRODIGAL(files_ch).faa

        MODELSEED(faa_ch.map { faa -> [faa.getBaseName(), faa] })

    } else {
        error "recon_with needs to be carveme, gapseq or modelseed (got: ${params.recon_with})"
    }

}
