#!/usr/bin/env nextflow

/* Nextflow script to run FAPROTAX

Usage: 

nextflow run modules/faprotax/main.nf -params-file modules/faprotax/faprotax.yaml 

*/

process FAPROTAX {

    tag "Annotate abundance table with FAPROTAX functional groups"

    publishDir "${params.outdir}/faprotax", mode: 'copy'
    container params.faprotax_version.toString() == "2" ? params.faprotax2_image : params.microbetag_image
    cpus params.threads

    input:
        path abundance_table
    
    output:
        path "functional_otu_table.tsv", emit: fapro_table
        path "sub_tables",               emit: fapro_subtables
        path "faprotax_report.txt",      optional: true
    
    script:
    if (params.faprotax_version.toString() == "2")
        // FAPROTAX v2 (faprotax.R): OTUs placed on the reference tree by their taxonomy. It reads a table by its
        // extension, so the table is linked with the one of its delimiter; its sub-tables are then written as
        // FAPROTAX 1.2.10 writes them (a comment line, and the abundance table's id column), for microbetag.
        """
        delim=\$(printf '%b' '${params.delimiter}')   # e.g. '\\t' as a tab
        ext=\$( [ "\$delim" = ',' ] && echo csv || echo tsv )
        ln -s ${abundance_table} input_table.\$ext
        id_col=\$(head -n 1 ${abundance_table} | cut -d "\$delim" -f 1)

        Rscript /opt/faprotax2/faprotax.R \
            -i input_table.\$ext \
            -d /opt/faprotax2/FAPROTAX_ref_data \
            -t "column:${params.taxonomy_colname}" \
            -k "\$id_col" \
            -o functional_otu_table.tsv \
            -b sub_tables \
            -r faprotax_report.txt \
            --Nthreads ${task.cpus} \
            -f

        for f in sub_tables/*; do
            sed -i "1s/^OTU/\$id_col/; 1i #" "\$f"
        done
        """
    else
        """
        python /opt/faprotax/collapse_table.py \
            -i ${abundance_table} \
            -o functional_otu_table.tsv \
            -g /opt/faprotax/FAPROTAX.txt \
            --table_delimiter '${params.delimiter}' \
            -c "#" \
            -d '${params.taxonomy_colname}' \
            -v \
            --force \
            -s sub_tables
        """
}


workflow{

    // Create a channel from input abundance table
    def abundance_ch = channel.fromPath("${params.abundance_table}")

    // Run process 
    FAPROTAX(abundance_ch)

}



