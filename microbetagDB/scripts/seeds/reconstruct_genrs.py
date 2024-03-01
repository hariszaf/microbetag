#!/usr/bin/env python3
import os
import re
import glob
import subprocess
import time
from modelseedpy import MSBuilder, MSGenome
import cobra

import multiprocessing

def split_list(input_list, chunk_size):
    """
    Split a list to sublists of a size.
    """
    return [input_list[i:i + chunk_size] for i in range(0, len(input_list), chunk_size)]


def handle_gtdb_id(gtdb_id):
    """
    Get PATRIC annotations for the GTDB genomes if available.
    The GTDB genomes can be found in the gtdb_ids.tsv file
    In the gtdb2patricIds.tsv file we keep the GTDB id along with its corresponding PATRIC id.
    Last, in the genomes_to_annotate.tsv file, we keep track of the GTDB ids for which there is no corresponding PATRIC annotation.
    """
    if gtdb_id in parsed_gtdb_ids:
        # print(gtdb_id, " parsed")
        return 1

    get_patic_genome_id_command = "".join(["p3-all-genomes --eq assembly_accession,", gtdb_id, " --attr genome_id --attr genome_name"])
    patric_genome_id_run = subprocess.check_output(get_patic_genome_id_command, shell=True, text=True)
    patric_genome_id = patric_genome_id_run.split("\n")[1].split("\t")[0]

    if patric_genome_id == "":

        if gtdb_id.startswith("GCF"):
            gtdb_id_gc_alt = gtdb_id.replace("GCF", "GCA")
        else:
            gtdb_id_gc_alt = gtdb_id.replace("GCA", "GCF")

        get_patic_genome_id_command = "".join(["p3-all-genomes --eq assembly_accession,", gtdb_id_gc_alt, " --attr genome_id --attr genome_name"])
        patric_genome_id_run = subprocess.check_output(get_patic_genome_id_command, shell=True, text=True)
        patric_genome_id = patric_genome_id_run.split("\n")[1].split("\t")[0]

    if patric_genome_id != "":
        gtdb_ids_to_patirc_ids = open("gtdb_genomes_reps_r207/gtdb2patricIds.tsv", "a")
        gtdb_ids_to_patirc_ids.write(gtdb_id + "\t" + patric_genome_id + "\n")
        print("new entry! ", gtdb_id)


    if patric_genome_id in patric_downs:
        print("already down ", gtdb_id)
        return 1

    wget_command = "".join(["wget -P patric_annotations/ -qN ftp://ftp.bvbrc.org/genomes/", patric_genome_id, "/", patric_genome_id, ".PATRIC.faa"])

    if os.system(wget_command) != 0:
        genomes_not_in_patric = open("gtdb_genomes_reps_r207/genomes_to_annotate.tsv", "a")
        genomes_not_in_patric.write(gtdb_id + "\n")
        print("new genome to annotate.")
    else:
        print("made it for gtdb id: ", gtdb_id)

    return 1


def rast_annotate_a_genome(filename):
    """
    Annotate a genome (.fna) file using the RAST server and the corresponding rast commands
    """
    start = time.time()

    # gunzip
    gzip_command = " ".join(["gunzip", filename])
    try:
        os.system(gzip_command)
    except:
        rast_annotate_a_genome(filename)

    decompressed_filename = filename[:-3]

    time.sleep(1)

    with open(decompressed_filename, 'r') as file:
        first_line = file.readline()


    # # Define the regex pattern to capture the desired text
    # pattern = r'(.*) (?:isolate|strain)'
    # match = re.search(pattern, first_line)
    # if match:
    #     result = match.group(1).strip()
    #     name = " ".join(result.split(" ")[1:])

    # pattern = r'(?:isolate|strain) (\S+)'
    # match = re.search(pattern, first_line)
    # if match:
    #     genome = match.group(1)

    name = " ".join(first_line.split(",")[0].split(" ")[1:])
    genome = filename.split("/")[-1].split("_")[2]

    gto_filename = "".join([decompressed_filename[:-4], ".gto"])

    # rast_create_genome_command
    rast_create_genome_command = " ".join(["rast-create-genome", "--scientific-name", name,
                            "--genetic-code", "11", "--domain", "Bacteria",
                            "--contigs", decompressed_filename,
                            "--genome-id", genome, ">", gto_filename])
    print("rast_create_genome_command:", rast_create_genome_command)
    os.system(rast_create_genome_command)

    # rast-process-genome
    gto_filename_2 = "".join([gto_filename, "_2"])
    rast_process_genome_command = " ".join(["rast-process-genome", "<", gto_filename, ">", gto_filename_2])
    print("\nrast_process_genome_command: ", rast_process_genome_command)
    try:
        os.system(rast_process_genome_command)
    except:
        pass

    # rast-export-genome protein_fasta
    faa_filename = "".join([gto_filename[:-4], ".faa"])
    rast_export_genome_command = " ".join(["rast-export-genome", "protein_fasta", "<", gto_filename_2, ">", faa_filename])
    print("\nrast_export_genome_command: ", rast_export_genome_command)
    try:
        os.system(rast_export_genome_command)
    except:
        pass

    # rast-export-genome seed_dir
    seed_dir = "".join([gto_filename[:-4], ".seed_dir.tar.gz"])
    rast_export_genome_dir_command = " ".join(["rast-export-genome", "seed_dir", "<", gto_filename_2, ">", seed_dir])
    print("\nrast_export_genome_dir_command:", rast_export_genome_dir_command)
    os.system(rast_export_genome_dir_command)

    end = time.time()

    print("RAST annotation was completed in", str((end-start)/60), "minutes.\nMove on with the ModelSEED part\n")


def modelseed_reconstruction(annotation_faa):
    """
    Build a Genome Scale reconstruction using ModelSEEDpy and the PATRIC annotations
    """
    # ModelSEED using the default b.f and complete medium, gapfill with the default algo
    msgenome = MSGenome.from_fasta(annotation_faa, split=' ')
    print("genome built: ", annotation_faa.split("/")[-1])
    with open(annotation_faa, 'r') as file:
        first_line = file.readline()
        model_id = first_line.split("[")[1]

        try:
            model_id = model_id.split("|")[0]
        except:
            pass

    print("about to reconstruct", annotation_faa.split("/")[-1])
    try:
        model = MSBuilder.build_metabolic_model(model_id = model_id,
                                            genome   = msgenome,
                                            index    = "0",
                                            classic_biomass = True,
                                            gapfill_model   = True,
                                            gapfill_media   = None,
                                            annotate_with_rast = True,
                                            allow_all_non_grp_reactions = True)

    except:
    #    modelseed_reconstruction(annotation_faa)
        return 0
    print("reconstruct ready")
    model_filename = "".join([annotation_faa[:-4], ".xml"])
    cobra.io.write_sbml_model(cobra_model = model, filename = model_filename)

    return 1



pwd = os.getcwd()

# patric_annotations = os.listdir("patric_annotations/")
# patric_annotations = [ "/".join([pwd, "patric_annotations", x]) for x in patric_annotations ]

# Load GTDB genomes ids and split them in lists of a thousand
gtdb_ids = open("gtdb_genomes_reps_r207/gtdb_ids.tsv", "r").readlines()
gtdb_ids = [line.strip() for line in gtdb_ids]
gtdb_ids_to_patirc_ids_per_thousand = split_list(gtdb_ids, 1000)

# PATRIC annotation files
patric_downs = os.listdir("patric_annotations")
patric_annotations_files = [ "/".join([pwd, "patric_annotations", x]) for x in patric_downs if x.endswith(".faa")]
patric_annotations_files = split_list(patric_annotations_files, 20)
patric_downs = [x.split(".PATRIC")[0] for x in patric_downs]

# 2-col file with GTDB id and its corresponding PATRIC id found
pairs = open("gtdb_genomes_reps_r207/gtdb2patricIds.tsv", "r").readlines()
parsed_gtdb_ids = [line.split("\t")[0] for line in pairs]




"""
# PART A: Run a pool to get all the PATRIC annotations

counter = 0
for a_thousand_ids in gtdb_ids_to_patirc_ids_per_thousand:

    number_processes = 10
    pool = multiprocessing.Pool(number_processes)
    missing_projects = pool.map(handle_gtdb_id, a_thousand_ids)
    pool.close()
    pool.join()
    counter += 1000
    print("We now have ", str(counter), " PATRIC annotations out of the ", str(len(gtdb_ids)))

"""



"""
# PART A.2: (out of this python script)

To get the genomes of the GTDB ids that were not annotated in the PATRIC database,
we ran the following 2 commands on the command line:

ncbi-genome-download -s genbank -A genomes_to_annotate.tsv --formats fasta -o genomes_not_in_patric -p 4 bacteria
ncbi-genome-download -s refseq -A genomes_to_annotate.tsv --formats fasta -o genomes_not_in_patric -p 4 bacteria

using the ncbi-genome-download: https://github.com/kblin/ncbi-genome-download
"""



# PART B: RAST annotate genomes that are not present in PATRIC database

#genome_ids = os.listdir("gtdb_genomes_reps_r207/genomes_not_in_patric")
#genome_ids = [ "/".join([pwd, "gtdb_genomes_reps_r207/genomes_not_in_patric", x]) for x in genome_ids if x.endswith(".fna.gz")]
#genome_ids_per_5 = split_list(genome_ids, 5)
#counter = 0



# for pentada in genome_ids_per_5:

#     number_processes = 10
#     pool = multiprocessing.Pool(number_processes)
#     missing_projects = pool.map(rast_annotate_a_genome, pentada)
#     pool.close()
#     pool.join()
#     counter += 5
#     print("We now have annotated ", str(counter), " genomes out of the ", str(len(genome_ids)))




# PART C: Reconstruct GEMs using modelSEEDpy
print(patric_annotations_files)
print("ready to go...")

for ahundred in patric_annotations_files:
    number_processes = 20
    pool = multiprocessing.Pool(number_processes)
    reconsrtuctions = pool.map(modelseed_reconstruction, ahundred)
    pool.close()
    pool.terminate()
    pool.join()
    #modelseed_reconstruction(acase)
    print("\n\nWe now have built a new reconstruction\n\n~~~~~~\n\n")




"""
# p3-all-genomes --eq assembly_accession,GCA_016713535.1 --attr genome_id --attr genome_name
# [NOT] curl -H 'Accept:application/protein+fasta' -H 'Content-Type:application/rqlquery+x-www-form-urlencoded' 'https://www.patricbrc.org/api/genome_feature/?eq(genome_id,83332.12)'

for i in $(ls *.xml) ; do mv ${i%.*}*  done/ ; done

"""
