//
// Subworkflow for seed complementarity pre-calculations for microbetag
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { readParamsFile; getInputFiles; sanitizeChannel; chunkFiles; SAFENAME_FILES; GUNZIP } from '../../modules/helpers'
include { GET_SEED_SETS; AGGREGATE_SEED_SETS; SCORES_COMPLS_PRECALC; AGGREGATE_SCORES_COMPLS } from '../../modules/seed_compl/seed_compl'
include { CARVE } from '../../modules/gem_recon/gem_recon'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow SEED_COMPLEMENTARITY {


    def has_gems = params.gems && file(params.gems).exists()
    gem_files_empty = true
    if (has_gems) {
        def gem_files_list = file(params.gems).list().findAll { 
             it.endsWith('.xml') || it.endsWith('.sbml') 
         }
        gem_files_empty = gem_files_list.isEmpty()
    }

    def species_ch
    def gems_ch

    if (gem_files_empty) {

        log.info "Genome-scale metabolic models were not provided and they will be reconstructed."

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

        // Reconstruct using software of user's choice
        if (params.recon_with =="carveme") {

            gems_ch = CARVE(infiles_chunk_ch, recon_in.is_faa)

        }

        species_ch = gems_ch.map { file -> file.baseName }

    } else {

        log.info "Genome-scale models were provided by the user and will be used. "

        // Channels for species GEMs and species names
        gems_ch    = Channel.fromPath("${params.gems}/*.xml")
        species_ch = Channel
            .fromPath("${params.gems}/*.xml")
            .map { file -> file.baseName }

    }

    def species_data_ch
    seedsets_avail = params.nonseeds_json && params.confidence_json && 
                        file(params.nonseeds_json).exists() && 
                        file(params.confidence_json).exists()

    if (!seedsets_avail) {

        log.info "Confidence and non-seed files were not provided and will be calculated"

        // Get seed and non-seed sets
        s = GET_SEED_SETS(gems_ch)

        // Aggregate seed and non-seed sets
        a = AGGREGATE_SEED_SETS(s.collect())

        // Combine scirpt and seed data for each species
        species_data_ch = species_ch
            .combine(a.nonseeds_json)
            .combine(a.confidence_json)

    } else {

        log.info "Previous calculated confidence and non-seed files were provided by the used and will be used."

        def nonseeds_json; def confidence_json
        nonseeds_json   = Channel.fromPath(params.nonseeds_json)
        confidence_json = Channel.fromPath(params.confidence_json)
        species_data_ch = species_ch
            .combine(nonseeds_json)
            .combine(confidence_json)
    }

    log.info "Calculate seed scores and complementarities using seed sets."

    // Calculate seed complementarity scores and extract complements
    e = SCORES_COMPLS_PRECALC(species_data_ch)

    // Aggregate compls and scores
    AGGREGATE_SCORES_COMPLS(e[0].collect(), e[1].collect())

}

