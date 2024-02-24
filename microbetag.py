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



# STEP1: phenotrex

suffixes = [".fa", ".fasta", ".gz"]
bin_files = get_files_with_suffixes(config.bins_path, suffixes)
bin_files_in_a_row = " ".join(bin_files)


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
    print(compute_genotype_command)
    print("Phenotrex genotype for the first time.")
    if os.system(compute_genotype_command) != 0:
        print("Try phenotrex genotype for the second time.")
        if os.system(compute_genotype_command) != 0:
            raise ValueError

folder_path = "microbetagDB/ref-dbs/phenDB/classes/"
phen_models = [os.path.join(folder_path, model) for model in os.listdir(folder_path)]



for model in phen_models:
    model_name =  os.path.basename(model)
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
        "".join([config.predictions_path, "/", model_name[:-4], ".prediction.tsv"])
    ]
    predict_trait_command = " ".join(predict_traits_params)
    print(predict_trait_command)

    if os.system(predict_trait_command) != 0:
        raise ValueError

