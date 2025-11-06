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



    // // 1️⃣ Create a channel of all .faa or .fa, .fasta files
    // def input_files = file("${params.input_files}")
    // def faa_files   = input_files.listFiles().findAll { it.name =~ /\.(faa|faa\.gz)$/ }
    // def nucl_files  = input_files.listFiles().findAll { it.name =~ /\.(fa|fna|fasta|fa\.gz|fna\.gz|fasta\.gz)$/ }

    // if (faa_files && !nucl_files) {
    //     log.info "Detected protein FASTA files (*.faa or *.faa.gz)"
    //     input_files_ch = sanitizeChannel("${params.input_files}/*.{faa,faa.gz}")
    //     is_faa = true
    // } else if (nucl_files && !faa_files) {
    //     log.info "Detected nucleotide FASTA files (*.fa, *.fna, *.fasta, etc.)"
    //     input_files_ch = sanitizeChannel("${params.input_files}/*.{fa,fasta,fna,fa.gz,fasta.gz,fna.gz}")
    //     is_faa = false
    // } else if (faa_files && nucl_files) {
    //     error "Mixed FASTA file types detected (both nucleotide and protein). Please separate them."
    // } else {
    //     error "No input FASTA files found in ${params.input_files}"
    // }

    // // sanitize filenames if for example filenames like BATCH:set1.fastq
    // prep = SAFENAME_FILES(input_files_ch)

    // // Decompress if .gz
    // decomp_files = GUNZIP(
    //     prep.map { orig_name, file ->
    //         def orig_name_decomp = orig_name.endsWith(".gz") ? orig_name[0..-4] : orig_name
    //         tuple(orig_name, orig_name_decomp, file)
    //     }
    // )

    // // 2️⃣ Collect all files into a list (for small/medium datasets)
    // input_list_ch = decomp_files.collect()

    // // 3️⃣ Compute chunk size and collate
    // infiles_chunk_ch = input_list_ch.flatMap { files ->

    //     log.info "Number of input files found: ${files.size()} "

    //     def chunk_size = Math.ceil(files.size() / params.max_forks) as int
    //     log.info "Chunk size = $chunk_size"

    //     // collate manually into sublists
    //     def chunks = []
    //     for (i = 0; i < files.size(); i += chunk_size) {
    //         chunks << files[i..Math.min(i+chunk_size-1, files.size()-1)]
    //     }

    //     return chunks
    // }