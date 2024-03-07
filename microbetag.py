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
# STEP1: phen annotations
# ----------------
print("\n PREDICTING PHENOTYPIC TRAITS \n")
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
print("\n PREDICTING ORFs \n")
for bin_fa in bin_files:
    bin_filename = os.path.basename(bin_fa)
    bin_id, extension = os.path.splitext(bin_filename)
    run_prodigal(bin_fa, bin_id, config.prodigal)

# ----------------
# STEP 3: KEGG annotation - using diting interface
# ----------------
print("\n KEGG ANNOTATION \n")
ko_list = os.path.join(config.kegg_db_dir, 'ko_list')
ko_dic = ko_list_parser(ko_list)

for bn in config.bin_filenames:
    bin_id, extension = os.path.splitext(bn)
    faa = os.path.join(config.prodigal, bin_id + '.faa')
    kegg_annotation(faa, bin_id, config.kegg_pieces_dir, config.kegg_db_dir, ko_dic, config.threads)

ko_merged_tab = os.path.join(config.kegg_annotations, 'ko_merged.txt')
bins_kos, pivot_df = merge_ko(config.kegg_pieces_dir, ko_merged_tab)

# ----------------
# STEP 4: Extract pathway complementarities
# ----------------
print("Exporting path complements....")
bin_kos_per_module, alt_to_gapfill, complements = export_pathway_complementarities(
    config,
    pivot_df
    )

alts_file = os.path.join(config.output_dir, "alts.json")
compl_file = os.path.join(config.output_dir, "pathCompls.json")
if not os.path.exists(alts_file):
    with open(alts_file, "w") as file:
        json.dump(alt_to_gapfill, file, cls=SetEncoder)
if not os.path.exists(compl_file):
    with open(compl_file, "w") as file:
        json.dump(complements, file, cls=SetEncoder)

# ----------------
# STEP 5: Build GENREs
# ----------------
print("\n BUILDING RECONSTRUCTIONS \n")
build_genres = build_genres(config)
build_genres.rast_annotate_genomes()
build_genres.modelseed_reconstructions()


# ----------------
# STEP 6: Phylomint
# ----------------
print("\n COMPUTING SEED SETS AND SCORES \n")
get_seed_sets(config)



# ----------------
# STEP 7: Export seed complementarities
# ----------------

seed_complements = export_seed_complementarities(config)
seed_complements.update()
# consider running again seed scores using update seed sets
# in this case, we should also edit the ConfidenceScore dictionary
# by removing seeds that were removed in the update()

seed_complements.module_related_seeds()

