#!/usr/bin/env python3

"""
Just like in the MGnify catalogues and the KEGG genomes cases, 
we need to have the KOs found per NCBI Taxonomy Id from the GTDB genomes, 
in a 2-column file, as shown here:

md:M00016       ko:K14267
md:M00034       ko:K00797
md:M00133       ko:K00797
"""

import os, sys

accessions      = []
accession_files = {}
mappings        = {}

metadata_file = open("GTDB_QUALITY_REPRESENTATIVE_GENOMES", "r")

for filename in os.listdir("../gtdb_genomes/ko_annotations"):
    accession = '_'.join([filename.split("_")[0], filename.split("_")[1]])
    accessions.append(accession)
    accession_files[accession] = filename

uniq_ncbi_ids = set()

for line in metadata_file: 

    line = line.split('\t')

    ncbi_genome_accession_id = line[54]
    ncbi_taxonomy_id         = line[77] 

    if ncbi_genome_accession_id in accessions:
        mappings[ncbi_genome_accession_id] = ncbi_taxonomy_id
        uniq_ncbi_ids.add(ncbi_taxonomy_id)

	# replace testingpy with ko_annotationsal 
        ncbi_path = os.path.join("../gtdb_genomes/testingpy/", ncbi_taxonomy_id)

        if os.path.isdir(ncbi_path) == False:
            os.mkdir(ncbi_path)

        faa_file        = "../gtdb_genomes/testingpy/" + accession_files[ncbi_genome_accession_id]
        kegg_annot_file = "../gtdb_genomes/testingpy/" + str(ncbi_taxonomy_id) + "/" + accession_files[ncbi_genome_accession_id]

        os.replace(faa_file, kegg_annot_file)

        kos_per_md = {}
        for line in open("../kegg_mappings/kegg_terms_per_module.tsv", "r"): 
            
            md = line.split("\t")[0]
            ko = line.split("\t")[1][:-1]
            if md not in kos_per_md:
                kos_per_md[md] = [ko[3:]]
            else: 
                kos_per_md[md].append(ko[3:])

        kos_per_mo_filename = kegg_annot_file.replace(accession_files[ncbi_genome_accession_id], ncbi_genome_accession_id + "_kos_related_to_mos.tsv")

        q = open(kos_per_mo_filename, "w")

        species_kos = set()
        for line in open(kegg_annot_file, "r"):

            line = line.split("\t")
            if len(line) == 1:
                continue
            else:

                ko_term =  line[-1][:-1]

                species_kos.add(ko_term)


        for ko_term in species_kos:

            for md, kos in kos_per_md.items():

                if ko_term in kos: 
                    q.write(md + "\t" + ko_term + "\n")

