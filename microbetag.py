"""
Aim:
    Perform the microbetag () approach annotating user's bins instead of mapping taxonomies against the
    microbetagDB representative genomes.

Input:
    - A folder with .fasta files of the corresponding bins
    - The abundance table of the bins across the samples
    - (optional) a co-occurrence network in a three-columns format

Output:
    -

Author:
    Haris Zafeiropoulos

"""
import os
import sys
import yaml
from config import Config
from utils import *

config_file = sys.argv[1]

with open(config_file, 'r') as yaml_file:
    config = Config(yaml.safe_load(yaml_file))

if config.bins_path is None:
    raise ValueError


# ----------------
# STEP1: phenotrex
# ----------------
suffixes = [".fa", ".fasta", ".gz"]
bin_files = get_files_with_suffixes(config.bins_path, suffixes)
bin_files_in_a_row = " ".join(bin_files)

# Build genotypes
compute_genotype_params = [ "phenotrex",
                            "compute-genotype",
                            "--out",
                            config.genotypes_file,
                            "--threads",
                            str(config.threads),
                            bin_files_in_a_row
                        ]
compute_genotype_command = " ".join(compute_genotype_params)

if not os.path.exists(config.genotypes_file):
    print("Phenotrex genotype for the first time.")
    if os.system(compute_genotype_command) != 0:
        print("Try phenotrex genotype for the second time.")
        if os.system(compute_genotype_command) != 0:
            raise ValueError
else:
    print("Genotypes already comptued.")

# Get predictions
folder_path = "microbetagDB/ref-dbs/phenDB/classes/"
phen_models = [os.path.join(folder_path, model) for model in os.listdir(folder_path)]

for model in phen_models:
    model_name =  os.path.basename(model)
    model_predictions_output = "".join([
        config.predictions_path, "/", model_name[:-4], ".prediction.tsv"
    ])
    if os.path.exists(model_predictions_output):
        print("Predictions already there for model:", model_name)
    else:
        predict_traits_params = [
            "phenotrex",
            "predict",
            "--genotype",
            config.genotypes_file,
            "--classifier",
            model,
            "--min_proba",
            str(config.min_proba),
            "--verb >",
            model_predictions_output
        ]
        predict_trait_command = " ".join(predict_traits_params)
        if os.system(predict_trait_command) != 0:
            raise ValueError

# ----------------
# STEP 2: prodigal - using diting interface
# ----------------

for bin_fa in bin_files:
    bin_filename = os.path.basename(bin_fa)
    bin_id, extension = os.path.splitext(bin_filename)
    run_prodigal(bin_fa, bin_id, config.prodigal)

# ----------------
# STEP 3: KEGG annotation - using diting interface
# ----------------

ko_list = os.path.join(config.kegg_db_dir, 'ko_list')
ko_dic = ko_list_parser(ko_list)

for bn in config.bin_filenames:
    bin_id, extension = os.path.splitext(bn)
    faa = os.path.join(config.prodigal, bin_id + '.faa')
    kegg_annotation(faa, bin_id, config.kegg_pieces_dir, config.kegg_db_dir, ko_dic, config.threads)

ko_merged_tab = os.path.join(config.kegg_annotations, 'ko_merged.txt')
merge_ko(config.kegg_pieces_dir, ko_merged_tab)


