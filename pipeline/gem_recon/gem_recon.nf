#!/usr/bin/env nextflow

/* Nextflow script to perform seed complementarity precalculations 

Usage: 

nextflow run seed_compl/seed_compl.nf -c seed_compl/seed_compl.config

or

nextflow run seed_compl/seed_compl.nf -params-file params/seed_compl.yaml
*/



include { readParamsFile } from '../helpers.nf'

// Only read default YAML if user didn't specify a params-file
if (!workflow.commandLine.contains('-params-file')) {
    println "[INFO] No params-file provided, loading default YAML..."
    def new_params = readParamsFile(params.paramsFile)
    params.putAll(new_params)
} else {
    println "[INFO] Using user-provided params-file, skipping default YAML."
}

println "Parameters after merge: ${params}"


params.max_forks = params.max_forks ?: 5


process carve {

    tag "GEM reconstruction with CarveMe"

    publishDir "${params.outdir}/reconstructions", mode: 'copy', overwrite: true
    container "carveme"
    containerOptions = "-v ${params.gurobi_lic}:/opt/gurobi/gurobi.lic"
    
    input:
    path in_chunk

    // Optional: uncomment to enable conditional execution
    when: 
    params.genre_reconstruction_with == "carveme"

    output:
    path "*.xml"

    script:
    """
    for f in ${in_chunk}; do

        if [ "${params.is_faa}" == "true" ]; then
            carve --solver gurobi --output "\$(basename \$f .faa).xml" \$f
        else
            base=\$(basename "$f") 
            base="\${base%.*}"
            carve --dna --solver gurobi --output "\${base}.xml" \$f
        fi

    done
    """
}

            // base="\${f%.fa}"
            // base="\${base%.fasta}"
            // base="\${base%.fna}"

process gapseq {

    tag "GEM reconstruction with gapseq"

    publishDir "${params.outdir}/reconstructions", mode: 'copy', overwrite: true
    container "gapseq"
    containerOptions = "-v ${params.cplex_lic}:/opt/cplex/cplex.lic"

    input:
    path in_chunk

    when:
    params.genre_reconstruction_with == "gapseq"

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

    // 1️⃣ Create a channel of all .faa or .fa, .fasta files
    def input_files_ch

    if (params.is_faa) {
        input_files_ch = Channel.fromPath("${params.input_files}/*.faa")
    } else {
        input_files_ch = Channel.fromPath("${params.input_files}/*.{fa,fasta,fna}")
    }

    // 1️⃣ Create a channel of all .faa files
    // input_files_ch = Channel.fromPath("${params.faa}/*.faa")

    // 2️⃣ Collect all files into a list (for small/medium datasets)
    input_list_ch = input_files_ch.collect()

    // 3️⃣ Compute chunk size and collate
    infiles_chunk_ch = input_list_ch.flatMap { files ->

        println "Number of input files found: ${files.size()} "

        def chunk_size = Math.ceil(files.size() / params.max_forks) as int
        println "Chunk size = $chunk_size"

        // collate manually into sublists
        def chunks = []
        for (i = 0; i < files.size(); i += chunk_size) {
            chunks << files[i..Math.min(i+chunk_size-1, files.size()-1)]
        }

        return chunks
    }

    infiles_chunk_ch.view()

    // 4️⃣ Pass chunks to your carve process
    carve(infiles_chunk_ch)

    // 4️⃣ Pass chunks to your gapseq process
    gapseq(infiles_chunk_ch)


}
