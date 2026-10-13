

// Mount CPLEX in the containers, if given (params.cplex); gapseq uses its default solver otherwise
def cplexBind() {
    if (!params.cplex) {
        return ""
    }
    return workflow.containerEngine == 'singularity'
        ? "-B ${params.cplex}:/opt/CPLEX"
        : "-v ${params.cplex}:/opt/CPLEX"
}


process GAPSEQ_FIND {

    tag "gapseq find"
    
    container "hariszaf/gapseq:1.4.0"
    publishDir path: { "${params.outdir}/reconstructions/${base}" }, mode: 'copy', overwrite: true

    input:
    tuple val(base), path(genome)

    output:
    tuple val(base), path("*-Reactions.tbl"), path("*-Pathways.tbl"), emit: find

    script:
    """
    base=\$(basename "${genome}")
    base="\${base%.*}"

    rxn_tbl=\${base}-${params.gapseq.find.p}-Reactions.tbl
    ptw_tbl=\${base}-${params.gapseq.find.p}-Pathways.tbl

    if [[ -f "\${rxn_tbl}" && -f "\${ptw_tbl}" ]]; then
        echo "[GAPSEQ FIND] Both -Reactions and -Pathways files were provided by the user. Skip gapseq find step." >> find.log
    else
        gapseq find \
            -p ${params.gapseq.find.p} \
            -t ${params.gapseq.find.t} \
            -b ${params.gapseq.find.b} \
            -i ${params.gapseq.find.i} \
            -c ${params.gapseq.find.c} \
            "${genome}"
        echo "[GAPSEQ FIND] Gapseq find completed for genome: \${base}." >> find.log
    fi 

    """
}


process GAPSEQ_FIND_TRANSPORTER {
    tag "gapseq find transporters"

    container "hariszaf/gapseq:1.4.0"
    publishDir path: { "${params.outdir}/reconstructions/${base}" }, mode: 'copy', overwrite: true

    input:
    tuple val(base), path(genome)

    output:
    tuple val(base), path("*-Transporter.tbl"), emit: trp_tbl
    path "transporter.log"

    script:
    """
    base=\$(basename "${genome}")
    base="\${base%.*}"

    echo -e "[GAPSEQ FIND-TRANSPORT] Run the find-transport script of gapseq -- searching for transporters... "  >> transporter.log

    trp_tbl=\${base}-Transporter.tbl

    if [[ -f "\${trp_tbl}" ]]; then 
        echo "[GAPSEQ FIND-TRANSPORT] Transporters file was provided by the user. Skip gapseq find-transporter step.." >> transporter.log

    else 
        gapseq find-transport \
            -b ${params.gapseq.find_transport.b} \
            -i ${params.gapseq.find_transport.i} \
            -c ${params.gapseq.find_transport.c} \
            "${genome}"
        echo "[GAPSEQ FIND-TRANSPORT] Gapseq find-transporters completed for genome: \${base}." >> transporter.log
    fi
    """
}


process GAPSEQ_DRAFT {

    tag "gapseq draft reconstruction"

    container "hariszaf/gapseq:1.4.0"
    containerOptions { cplexBind() }

    publishDir path: { "${params.outdir}/reconstructions/${base}" }, mode: 'copy', overwrite: true
    

    input:
    tuple val(base),
        path(genome),
        path(rxn_tbl), // "*-Reactions.tbl"
        path(pwt_tbl), 
        path(trp_tbl)

    output:
    tuple val(base), 
        path("*-draft.RDS"), 
        path("*-rxnWeights.RDS"), 
        path("*-rxnXgenes.RDS"), 
        emit: draft_recon_files
    path("*-draft.xml"), emit: draft_xml

    script:
    """
    # ---------
    # Check if CPLEX is available, and if so, install cobrarCPLEX
    # ---------

    if [ -d /opt/CPLEX ] && [ "\$(ls -A /opt/CPLEX)" ]; then
        echo ">>>>> /opt/CPLEX is mounted and not empty. Let us nstall cobrarCPLEX .. " >> emm
        unzip /opt/Rsrc/cobrarCPLEX.zip -d /opt/Rsrc/
        R CMD INSTALL /opt/Rsrc/cobrarCPLEX-main --configure-args="--with-cplex-dir=/opt/CPLEX/cplex"
    else
        echo ">>>>> /opt/CPLEX is missing or empty, skipping cobrarCPLEX installation" >> emm
    fi


    echo "[GAPSEQ DRAFT]: Starting draft reconstruction" > draft.log

    if [[ ${params.gapseq.draft.b} == "" || ${params.gapseq.draft.b} == "auto" ]]; then 
        gapseq draft \
            -r ${rxn_tbl} \
            -p ${pwt_tbl} \
            -t ${trp_tbl} \
            -c "${genome}" \
            -u ${params.gapseq.draft.u} \
            -l ${params.gapseq.draft.l}
    else
        gapseq draft \
            -r ${rxn_tbl} \
            -p ${pwt_tbl} \
            -t ${trp_tbl} \
            -b ${params.gapseq.draft.b} \
            -u ${params.gapseq.draft.u} \
            -l ${params.gapseq.draft.l}
    fi
    """
}


process GAPSEQ_FILL {

    tag "gapseq gapfill"
    
    container "hariszaf/gapseq:1.4.0"
    containerOptions { cplexBind() }

    publishDir path: { "${params.outdir}/reconstructions/${base}" }, mode: 'copy', overwrite: true

    input:
    tuple val(base),
        path(draft_rds),
        path(rxn_wghts),
        path(rxn_genes),
        path(medium)

    output:
    path "${base}.RDS", emit: filled_rds
    path "${base}.xml", emit: filled_xml

    script:
    """
    # ---------
    # Check if CPLEX is available, and if so, install cobrarCPLEX
    # ---------

    if [ -d /opt/CPLEX ] && [ "\$(ls -A /opt/CPLEX)" ]; then
        echo ">>>>> /opt/CPLEX is mounted and not empty. Let us nstall cobrarCPLEX .. " >> emm
        unzip /opt/Rsrc/cobrarCPLEX.zip -d /opt/Rsrc/
        R CMD INSTALL /opt/Rsrc/cobrarCPLEX-main --configure-args="--with-cplex-dir=/opt/CPLEX/cplex"
    else
        echo ">>>>> /opt/CPLEX is missing or empty, skipping cobrarCPLEX installation" >> emm
    fi

    echo "[GAPSEQ FILL]: Starting gapfilling" > fill.log

    gapseq fill \
        -m ${draft_rds} \
        -c ${rxn_wghts} \
        -g ${rxn_genes} \
        -n ${medium} \
        -b ${params.gapseq.medium.b}

    """
}


// A file of a previous gapseq run (params.gapseq.prev_tbl), e.g. gapseqPrev('bin_101', '-draft.RDS'): in the folder
// itself, or in a folder per genome, as this module publishes them (reconstructions/bin_101/bin_101-draft.RDS)
def gapseqPrev(base, suffix) {
    def prev = params.gapseq?.prev_tbl
    if (!prev) {
        return null
    }
    def flat = file("${prev}/${base}${suffix}")
    return flat.exists() ? flat : file("${prev}/${base}/${base}${suffix}")
}

def hasGapseqPrev(base, suffixes) {
    return suffixes.every { suffix -> def f = gapseqPrev(base, suffix); f != null && f.exists() }
}


workflow GAPSEQ {

    take:
    // [base, genome] of each genome
    genomes_with_base

    main:
    def p = params.gapseq.find.p
    def find_files  = ["-${p}-Reactions.tbl", "-${p}-Pathways.tbl"]
    def trp_files   = ["-Transporter.tbl"]
    def draft_files = ["-draft.RDS", "-rxnWeights.RDS", "-rxnXgenes.RDS"]

    // gapseq find (reactions and pathways), unless available from a previous run or not needed (draft available)
    def find_prev_ch = genomes_with_base
        .filter { base, _genome -> hasGapseqPrev(base, find_files) }
        .map { base, _genome -> [base, gapseqPrev(base, find_files[0]), gapseqPrev(base, find_files[1])] }
    def find_ch = find_prev_ch.concat(
        GAPSEQ_FIND(genomes_with_base.filter { base, _genome -> !hasGapseqPrev(base, find_files) && !hasGapseqPrev(base, draft_files) }).find
    )

    // gapseq find-transport, likewise
    def trp_prev_ch = genomes_with_base
        .filter { base, _genome -> hasGapseqPrev(base, trp_files) }
        .map { base, _genome -> [base, gapseqPrev(base, trp_files[0])] }
    def trp_ch = trp_prev_ch.concat(
        GAPSEQ_FIND_TRANSPORTER(genomes_with_base.filter { base, _genome -> !hasGapseqPrev(base, trp_files) && !hasGapseqPrev(base, draft_files) }).trp_tbl
    )

    // gapseq draft, unless available from a previous run: [base, draft.RDS, rxnWeights.RDS, rxnXgenes.RDS]
    def draft_prev_ch = genomes_with_base
        .filter { base, _genome -> hasGapseqPrev(base, draft_files) }
        .map { base, _genome -> [base] + draft_files.collect { suffix -> gapseqPrev(base, suffix) } }
    def draft_input = genomes_with_base
        .filter { base, _genome -> !hasGapseqPrev(base, draft_files) }
        .join(find_ch)
        .join(trp_ch)
    def draft = GAPSEQ_DRAFT(draft_input)
    def draft_ch = draft_prev_ch.concat(draft.draft_recon_files)

    // gapseq fill, on the medium file; without one, the models are the drafts
    def med_file = params.gapseq?.medium?.file ?: null
    def models
    if (med_file) {
        def med_ch = channel.value(file(med_file, checkIfExists: true))
        models = GAPSEQ_FILL(draft_ch.combine(med_ch)).filled_xml
    } else {
        log.warn "No medium file for gapseq (gapseq.medium.file): the models are not gap-filled"
        def draft_xml_prev = genomes_with_base
            .filter { base, _genome -> hasGapseqPrev(base, draft_files + ["-draft.xml"]) }
            .map { base, _genome -> gapseqPrev(base, "-draft.xml") }
        models = draft_xml_prev.concat(draft.draft_xml)
    }

    emit:
    model = models
}


workflow {

    def genomes = channel.fromPath("${params.recon_files}/*.{fa,fna,fasta}")
    GAPSEQ(genomes.map { genome -> [genome.getBaseName(), genome] })
}
