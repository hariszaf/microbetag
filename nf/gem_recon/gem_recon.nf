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

    tag "GEM reconstruction with CarveMe using .faa files as input"

    publishDir "${params.outdir}/reconstructions", mode: 'copy', overwrite: true
    container "carveme"
    containerOptions = "-v ${params.gurobi_lic}:/opt/gurobi/gurobi.lic"
    
    input:
    path faa_chunk

    // Optional: uncomment to enable conditional execution
    when: 
    params.genre_reconstruction_with == "carveme" && params.faa != null

    output:
    path "*.xml"

    script:
    """
    for f in ${faa_chunk}; do
        carve --solver gurobi --output "\$(basename \$f .faa).xml" \$f
    done
    """
}
    // export GRB_LICENSE_FILE=/opt/gurobi/gurobi.lic



workflow {

    genomes_ch   = Channel.fromPath("${params.genomes}")

    // 1️⃣ Create a channel of all .faa files
    faa_ch = Channel.fromPath("${params.faa}/*.faa")

    // 2️⃣ Collect all files into a list (for small/medium datasets)
    faa_list_ch = faa_ch.collect()

    // 3️⃣ Compute chunk size and collate
    faa_chunk_ch = faa_list_ch.flatMap { files ->

        println "Found ${files.size()} .faa files"

        def chunk_size = Math.ceil(files.size() / 5.0) as int
        println "Chunk size = $chunk_size"

        // collate manually into sublists
        def chunks = []
        for (i = 0; i < files.size(); i += chunk_size) {
            chunks << files[i..Math.min(i+chunk_size-1, files.size()-1)]
        }

        return chunks
    }

    // 4️⃣ Debug
    faa_chunk_ch.view { println "Chunk: $it" }

    // 5️⃣ Pass chunks to your carve process
    carve(faa_chunk_ch)


}




