#!/usr/bin/env nextflow

/* Nextflow script to perform Pathway complementarity precalculations 

Usage: 

nextflow run pathway_compl/pathway_compl.nf -c pathway_compl/kofam.config

or

nextflow run pathway_compl/pathway_compl.nf -params-file params/pathway_compl.yaml
*/

include { READPARAMSFILE } from '../helpers.nf'

// Only read default YAML if user didn't specify a params-file
if (!workflow.commandLine.contains('-params-file')) {
    println "[INFO] No params-file provided, loading default YAML..."
    def new_params = READPARAMSFILE(params.paramsFile)
    params.putAll(new_params)
} else {
    println "[INFO] Using user-provided params-file, skipping default YAML."
}

println "Parameters after merge: ${params}"


params.ko_output_file = params.ko_output_file ?: null

if( !params.ko_output_file ) {
    exit 1, "❌ The parameter 'ko_output_file' is required. Example: nextflow run main.nf --ko_output_file myfile.txt"
}



process PC_PRECALC {
    tag "Pathway complementarity precalculations"

    publishDir "${params.outdir}/pathway_compl", mode: 'copy'
    container "hariszaf/microbetag-nf:0.1.0"

    input:
    path ko_merged_ch
    path path_compl_sc_ch

    output:
    path "${params.alts_file}", emit: alts
    path "${params.pc_file}", emit: pcompls

    script:
    """
    python ${path_compl_sc_ch} \
        ${params.alts_file} \
        ${params.pc_file} \
        ${ko_merged_ch} \
        ${params.tinyurl} \
        ${params.threads}
    """
}
// ${params.pc_percent} \

process PC_EXTEND {
    tag ""

    publishDir "${params.outdir}/pathway_compl", mode: 'copy'
    container "hariszaf/microbetag-nf:0.1.0"

    input:
    path extend_sc
    path pcompls

    output:
    path "${pcompls.baseName}_ext.json", emit: pcompls_ext
    
    script:
    """
    python ${extend_sc} \
        ${pcompls} \
        ${params.pc_percent} \
        ${params.threads}

    mv pathway_complements_extended.json ${pcompls.baseName}_ext.json
    """

}


workflow {

    def ko_merged_ch     = Channel.fromPath(params.ko_output_file)
    def path_compl_sc_ch = Channel.fromPath("modules/pathway_compl/pathway_compl.py")
    def extend_sc_ch     = Channel.fromPath("modules/pathway_compl/extend.py")

    // (alts, pcompls) = PC_PRECALC(ko_merged_ch, path_compl_sc_ch)
    PC_PRECALC(ko_merged_ch, path_compl_sc_ch)

    // PC_EXTEND(alts, pcompls)
    PC_EXTEND(extend_sc_ch, PC_PRECALC.out.pcompls)

}
