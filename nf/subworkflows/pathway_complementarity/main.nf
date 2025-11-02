//
// Subworkflow for pathway complementarity pre-calculations for microbetag
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { READPARAMSFILE } from './modules/helpers.nf'
include { PRODIGAL } from './modules/prodigal/prodigal.nf'
include { HMMSEARCH; MERGE_HMMOUT; KO_ANNOTATE } from './modules/kofam/kofam.nf'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PATHWAY_COMPLEMENTARITY{

    take:
    genomes
    faa
    ko_merged

    main:

    // Input (annotated) genome files 
    def genomes_ch = Channel.fromPath("${params.genomes}/*")

    // KOFAM db
    def ko_list_ch      = Channel.fromPath("${params.kegg_list}")
    def hmm_prof_ch     = Channel.fromPath("${params.hmm_profiles}")


    // Create a channel from input genomes
    def hmmsearch_sc_ch = Channel.fromPath('modules/kofam/kofam.sh')


    def faa_ch


    // PATHWAY COMPLEMENTARITY STEP
    def alts_file
    def pc_file
    def ko_merged_ch = Channel.fromPath(params.ko_output_file, checkIfExists: true)


    // Check if we can skip PC_PRECALC
    def skip_precalc = params.alts_file && params.pc_file && 
                      file(params.alts_file).exists() && 
                      file(params.pc_file).exists()
    
    // get_ko2contig should only run if we're NOT skipping precalc
    def get_ko2contig = !skip_precalc && params.ko_output_file && file(params.ko_output_file).exists()


   // Check if faa

    // Check if faa_dir exists and contains files
    def faa_dir = params.faa_dir ? file(params.faa_dir) : null
    




    // path_compl.nf module
    if (skip_precalc) {
        alts_file = Channel.fromPath(params.alts_file, checkIfExists: true)
        pc_file   = Channel.fromPath(params.pc_file, checkIfExists: true)
        log.info "✓ Using existing files, skipping PC_PRECALC"

    } else {


        if (get_ko2contig) {

                // faa_ch = PRODIGAL(genomes_ch)  // returns a channel of *.faa files


            // Combine the single channels with all genomes
            def inputs_ch = faa_ch
                .combine(hmmsearch_sc_ch)
                .combine(ko_list_ch)
                .combine(hmm_prof_ch)

            // Run HMMSEARCH for each genome/bin
            hmmout_ch = HMMSEARCH(inputs_ch)

            // Collect all emitted hmmout dirs (waits for all tasks to finish)
            merged_input_ch = hmmout_ch.collect()

            // Merge all hmmout results into a single file
            MERGE_HMMOUT(merged_input_ch, merge_sc_ch)

        }


        PC_PRECALC(ko_merged_ch)
        alts_file = PC_PRECALC.out.alts
        pc_file   = PC_PRECALC.out.pcompls
        log.info "○ Running PC_PRECALC to generate files"
    }
    
    // Continue pipeline
    PC_EXTEND(pc_file)



}




    // // Check if faa_dir exists and contains files
    // def faa_dir = params.faa_dir ? file(params.faa_dir) : null

    // def faa_ch

    // if (faa_dir?.exists() && faa_dir?.isDirectory() && faa_dir?.listFiles()?.find { it.name.endsWith('.faa') }) {
    //     println "Using existing FAA files from ${params.faa_dir}"
    //     faa_ch = Channel.fromPath("${params.faa_dir}/*.faa")
    // } else {
    //     println "FAA files missing or empty, running ORFS/PRODIGAL"
    //     def genomes_ch = Channel.fromPath("${params.genomes}/*")
    //     faa_ch = PRODIGAL(genomes_ch)  // returns a channel of *.faa files
    // }

