#!/usr/bin/env nextflow

/* Nextflow script to perform seed complementarity analysis on GEMs

Usage: 

nextflow run seed_compl/seed_compl.nf -c seed_compl/seed_compl.config
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


process get_seed_sets {

    tag "Calculate seed and non-seed sets based on GEM reconstructions."

    publishDir "${params.outdir}/seed_compl/sets", mode: 'copy'
    container "hariszaf/microbetag-nf:0.1.0"

    input:
    tuple path(sets_sc), path(gem)

    output:
    path "${gem.baseName}_sets.json", emit: seed_sets_json

    script:
    """    
    python ${sets_sc} ${gem} ${params.namespace} "${gem.baseName}_sets.json"
    """

}

process aggregate_seed_sets {

    tag "Aggregate seed and non-seed sets, as JSON files, and filter them based on KEGG modules, pickle files."

    publishDir "${params.outdir}/seed_compl", mode: 'copy'
    container "hariszaf/microbetag-nf:0.1.0"

    input:
    path seed_sets_json
    path build_pkls_sc

    output:
    path "seeds.json", emit: seeds_json
    path "nonseeds.json", emit: nonseeds_json
    path "confidence.json", emit: confidence_json
    path "keggm_seeds.pkl", emit: module_related_seeds
    path "keggm_nonseeds.pkl", emit: module_related_nonseeds

    script:
    """
    jq -s '.' ${seed_sets_json.join(' ')} > all_seed_sets.json
    jq 'map({ (.base_name): .seed_set }) | add' all_seed_sets.json > seeds.json
    jq 'map({ (.base_name): .non_seed_set }) | add' all_seed_sets.json > nonseeds.json
    jq 'map({ (.base_name): .seed_set_confidence }) | add' all_seed_sets.json > confidence.json

    python ${build_pkls_sc}
    """
}

process scores_and_compl_precalc {

    tag "Calculate seed complementarity scores and extract complements."

    publishDir "${params.outdir}/seed_compl/per_species", mode: 'copy'
    container "hariszaf/microbetag-nf:0.1.0"

    input:
    tuple val(species), path(extract_sc), path(nonseeds_json), path(confidence_json)

    output:
    path "${species}_scores.tsv", emit: sp_seed_scores_json
    path "${species}_compls.json", emit: sp_seed_compl_json

    script:
    """
    python ${extract_sc} \
        "${species}" \
        "${nonseeds_json}" \
        "${confidence_json}" \
        "${params.kegg_modules_only}" \
        "${species}_scores.tsv" \
        "${species}_compls.json"
    """
}


process aggregate_scores_compls {
    
    tag "Aggregate the per species seed scores and complements to global files"

    publishDir "${params.outdir}/seed_compl", mode: 'copy'
    container "hariszaf/microbetag-nf:0.1.0"

    input:
    path sp_seed_scores_json
    path sp_seed_compl_json

    output:
    path compls_js_outfile,  emit: tmp_compls
    path compls_pkl_outfile, emit: seed_compls
    path scores_outfile,     emit: seed_scores

    script:
    compls_js_outfile  = "all_compls.json"
    compls_pkl_outfile = "seed_compls.pkl"
    scores_outfile     = "seed_scores.tsv"
    """ 
    jq -n '
      reduce inputs as \$f ({}; . + {(\$f|input_filename|capture("(?<key>[^/]+)_compls\\\\.json\$").key): \$f})
    ' *_compls.json > ${compls_js_outfile}

    echo -e 'nodeA\\tnodeB\\tCompetitionScore\\tCooperationScore' > header
    cat *.tsv >> scores
    cat header scores > ${scores_outfile}
    
    python - <<'EOF'
    import json
    import pickle
    import pandas as pd

    compls_js_outfile  = "${compls_js_outfile}"
    compls_pkl_outfile = "${compls_pkl_outfile}"

    with open(compls_js_outfile) as f:
        compls_dict = json.load(f)

    df = pd.DataFrame.from_dict(compls_dict)

    # Identify only the float columns
    float_cols = df.select_dtypes(include="float").columns

    # Replace NaNs with [] only in those columns
    df[float_cols] = df[float_cols].where(df[float_cols].notna(), [[]])

    # Attention! We need to get df.T. Otherwise we get the source as target and the other way around !
    with open(compls_pkl_outfile, "wb") as f:
        pickle.dump(df.T, f)
    EOF
    """
}


workflow {

    // Channels for species GEMs and species names
    gems_ch    = Channel.fromPath("${params.gems}/*.xml")
    species_ch = Channel
        .fromPath("${params.gems}/*.xml")
        .map { file -> file.baseName }

    // Channels for seed complementarity scripts
    sets_sc_ch    = Channel.fromPath("modules/seed_compl/get_seed_sets.py")
    extract_sc_ch = Channel.fromPath("modules/seed_compl/seed_compls.py")
    pkl_sc_ch     = Channel.fromPath("modules/seed_compl/build_pkls.py")

    // Get seed and non-seed sets
    s = get_seed_sets(sets_sc_ch.combine(gems_ch))

    // Aggregate seed and non-seed sets
    a = aggregate_seed_sets(s.collect(), pkl_sc_ch)

    // Combine scirpt and seed data for each species
    species_data_ch = species_ch
        .combine(extract_sc_ch)
        .combine(a.nonseeds_json)
        .combine(a.confidence_json)

    // Calculate seed complementarity scores and extract complements
    e = scores_and_compl_precalc(species_data_ch)

    // Aggregate compls and scores
    aggregate_scores_compls(e[0].collect(), e[1].collect())
}

