#!/usr/bin/env nextflow

/* Taxonomy assignment of the sequences of an amplicon abundance table with IDTAXA (DECIPHER),
against the 16S rRNA gene sequences of the GTDB genomes.

The abundance table has the sequence ids in its first column, the samples, and the sequences in its
last column. The output has the GTDB taxonomy in place of the sequences, as microbetag expects it.

Usage:

nextflow run modules/taxonomy/main.nf -params-file modules/taxonomy/taxonomy.yaml
*/

process TO_FASTA {

    tag "Write the sequences of the abundance table as FASTA"

    container params.microbetag_image

    input:
        path abundance_table

    output:
        path "sequences.fasta", emit: fasta

    script:
        """
        prep_taxonomy.py fasta ${abundance_table} sequences.fasta
        """
}


process IDTAXA {

    tag "Classify the sequences with IDTAXA against the GTDB 16S training set"

    container "quay.io/biocontainers/bioconductor-decipher:3.2.0--r44h3df3fcb_1"
    cpus params.threads ?: 2

    input:
        path fasta
        path training_set
        path classify_script

    output:
        path "assignments.tsv", emit: assignments

    script:
        """
        Rscript ${classify_script} ${fasta} ${training_set} ${params.idtaxa_threshold ?: 60} assignments.tsv
        """
}


process MERGE_TAXONOMY {

    tag "Replace the sequences of the abundance table with their GTDB taxonomy"

    publishDir "${params.outdir}", mode: 'copy'
    container params.microbetag_image

    input:
        path abundance_table
        path assignments

    output:
        path "GTDB_tax_assigned_abundance_table.tsv", emit: table

    script:
        """
        prep_taxonomy.py merge ${abundance_table} ${assignments} GTDB_tax_assigned_abundance_table.tsv
        """
}


workflow TAXONOMY {

    take:
        abundance_table
        training_set

    main:
        classify_script = channel.fromPath("${moduleDir}/../../scripts/taxonomy/classify.R", checkIfExists: true)

        fasta       = TO_FASTA(abundance_table)
        assignments = IDTAXA(fasta.fasta, training_set, classify_script)
        merged      = MERGE_TAXONOMY(abundance_table, assignments.assignments)

    emit:
        table = merged.table
}


workflow {

    TAXONOMY(
        channel.fromPath(params.abundance_file, checkIfExists: true),
        channel.fromPath(params.gtdb_16s_training_set, checkIfExists: true)
    )
}
