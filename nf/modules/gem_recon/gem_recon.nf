#!/usr/bin/env nextflow

/* Nextflow script to perform seed complementarity precalculations 

Usage: 

nextflow run seed_compl/seed_compl.nf -c seed_compl/seed_compl.config

or

nextflow run seed_compl/seed_compl.nf -params-file params/seed_compl.yaml
*/

include  { getInputFiles; sanitizeChannel; chunkFiles; SAFENAME_FILES; GUNZIP } from '../helpers.nf'


process CARVE {

    tag "GEM reconstruction with CarveMe"

    publishDir "${params.outdir}/reconstructions", mode: 'copy', overwrite: true
    container "hariszaf/carveme:1.6.6"
    containerOptions = "-v ${params.gurobi_lic}:/opt/gurobi/gurobi.lic"
    
    input:
        path in_chunk
        val is_faa

    // Optional: uncomment to enable conditional execution
    when:
        params.recon_with == "carveme"

    output:
        path "*.xml"

    script:
        """
        for f in ${in_chunk}; do

            base=\$(basename "\$f") 
            base="\${base%.*}"

            if [ ${is_faa} == "true" ]; then
                carve --solver gurobi --output "\${base}.xml" \$f
            else
                carve --dna --solver gurobi --output "\${base}.xml" \$f
            fi

        done
        """
}


process GAPSEQ {

    tag "GEM reconstruction with gapseq"

    publishDir "${params.outdir}/reconstructions", mode: 'copy', overwrite: true
    container "gapseq"
    containerOptions = "-v ${params.cplex_lic}:/opt/cplex/cplex.lic"

    input:
        path in_chunk

    when:
        params.recon_with == "gapseq"

    output:
        path ".xml"

    script:
        """
        for f in ${in_chunk}; do 

            gapseq doall \$f

        done
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

    // Use the COLLECTED channel for chunking, not the original decomp channel
    infiles_chunk_ch = chunkFiles(decomp_ch, params.max_forks)

    // Pass chunks to your carve process
    CARVE(infiles_chunk_ch, recon_in.is_faa)

    // // Pass chunks to your gapseq process
    // GAPSEQ(infiles_chunk_ch)

}
