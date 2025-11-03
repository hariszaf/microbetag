//
// Subworkflow for pathway complementarity pre-calculations for microbetag
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { READPARAMSFILE } from '../../modules/helpers'
include { PRODIGAL } from '../../modules/prodigal/prodigal'
include { HMMSEARCH; MERGE_HMMOUT } from '../../modules/kofam/kofam'
include { PC_PRECALC; PC_EXTEND } from '../../modules/pathway_compl/pathway_compl'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/



workflow PATHWAY_COMPLEMENTARITY {

    // PATHWAY COMPLEMENTARITY STEP
    def ko2contig_ch
    def alts_file
    def pc_file

    // Extract ko2contig
    def faa_ch
    def ko_list_ch
    def hmm_prof_ch
    def hmmsearch_sc_ch
    def hmm_in_ch
    
    def has_faa_dir

    // Genome annotation 
    def genomes_ch
    def prodigal_output

    def prec_res


    // Check if we can skip ALL upstream processing
    def skip_all_upstream = params.alts_file && params.pc_file && 
                           file(params.alts_file).exists() && 
                           file(params.pc_file).exists()

    if (skip_all_upstream) {

        // Skip everything - just use the provided files
        alts_file = file(params.alts_file)
        pc_file   = file(params.pc_file)
        log.info "✓ Using existing alts_file and pc_file, skipping ALL upstream processing"

        // Continue pipeline with the generated/already available files
        PC_EXTEND(pc_file)


    } else {

        // Need to run some or all of the upstream pipeline
        log.info "○ Running upstream pipeline to generate alts_file and pc_file"

        def has_ko_output = params.ko_merged && file(params.ko_merged).exists()

        // If ko to conting exists, you can run the generate alts_file and pc_file directly 
        if (has_ko_output) {

            // Use existing KO output file, skip annotation
            ko2contig_ch = Channel.fromPath(params.ko_merged, checkIfExists: true)
            log.info "✓ Using existing KO output file: ${params.ko_merged}"
        
        } else {

            log.info "○ Building the KO to contig file"

            // KOFAM db - validate required inputs
            if (!params.kegg_list || !file(params.kegg_list).exists()) {
                throw new Exception("KEGG list file does not exist: ${params.kegg_list}")
            }
            if (!params.hmm_profiles || !file(params.hmm_profiles).exists()) {
                throw new Exception("HMM profiles file does not exist: ${params.hmm_profiles}")
            }
            
            // Define KOFAM channels only when needed
            ko_list_ch      = Channel.fromPath(params.kegg_list, checkIfExists: true)
            hmm_prof_ch     = Channel.fromPath(params.hmm_profiles, checkIfExists: true)
            hmmsearch_sc_ch = Channel.fromPath('modules/kofam/kofam.sh')

            // Handle FAA files
            has_faa_dir = params.faa_dir && file(params.faa_dir).exists()
            faa_files_empty = true
            if (has_faa_dir) {
                def faa_files_list = file(params.faa_dir).list().findAll { it.endsWith('.faa') }
                faa_files_empty = faa_files_list.isEmpty()
            }

            if (has_faa_dir && !faa_files_empty) {

                // Use existing FAA files
                log.info "✓ Using existing FAA files from: ${params.faa_dir}"
                faa_ch = Channel.fromPath("${params.faa_dir}/*.faa", checkIfExists: true)
            
            } else {

                log.info "○ Running genome annotation to get ORFs with Prodigal "

                genomes_ch = Channel.fromPath("${params.genomes}/*.{fa, fasta}", checkIfExists: true)
                prodigal_output = PRODIGAL(genomes_ch)
                faa_ch =  prodigal_output.faa

            }

            // Combine the single channels with all faa_dir
            hmm_in_ch = faa_ch
                .combine(hmmsearch_sc_ch)
                .combine(ko_list_ch)
                .combine(hmm_prof_ch)
            
            // Run HMMSEARCH for each genome/bin
            def hmmout_ch = HMMSEARCH(hmm_in_ch)

            // Collect all emitted hmmout dirs (waits for all tasks to finish)
            def merged_input_ch = hmmout_ch.collect()

            // Merge all hmmout results into a single file
            def merge_res = MERGE_HMMOUT(merged_input_ch)
            ko2contig_ch = merge_res.ko2contig            
        }

        // Run PC_PRECALC with the KO data (either existing or generated)
        def pc_precalc_result = PC_PRECALC(ko2contig_ch)

        PC_EXTEND(pc_precalc_result.pcompls)

    }
}



