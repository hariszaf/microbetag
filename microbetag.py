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
    config = Config(yaml.safe_load(yaml_file), config_file)

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
    print("Genotypes already computed.")

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
if not os.path.exists(config.phylomint_scores):
    run_phylomint(config)
else:
    print("Seed scores already computed.")

# ----------------
# STEP 7: Export seed complementarities
# ----------------
print("\n EXORTING SEED COMPLEMENTS \n")
seed_complements = export_seed_complementarities(config)

"""
[NOTE]:consider running again "seed scores" (PhyloMint) using update seed sets
in this case, we should also edit the ConfidenceScore dictionary
by removing seeds that were removed in the update()
"""
if not os.path.exists(seed_complements.updated_seed_sets):
    seed_complements.update()
else:
    print("Seed sets already updated.")

if not os.path.exists(seed_complements.module_seeds):
    seed_complements.module_related_seeds()
else:
    print("Seed and non seed sets with compounds related to KEGG modules already retrieved.")

if not os.path.exists(seed_complements.seed_complements):
    seed_complements.export_seed_complements()
else:
    print("Seed complements already exported.")


# ----------------
# STEP 8: Build network if not available
# ----------------

if config.network is None:
    ensure_flashweave_format(conf=config)
    flashweave_params = [
        "julia",
        config.flashweave_script,
        config.mount,
        config.flashweave_abd_table,
        config.sensitive,
        config.heterogeneous,
        config.metadata,
        config.metadata_file
    ]
    flashweave_command = " ".join(flashweave_params)

    print("Run FlashWeave")
    print(flashweave_command)
    if os.system(flashweave_command) != 0:
        e = """ \
            FlashWeave failed.
            Please check on the content of your abundance table and the format of your metadata file if one provided.
            We suggest you use FlashWeave or any other software to get the network and then perform microbetag providing the network returned.
            You may find an implementation of FlashWeave in a Docker image in microbetag's preprocessing image:
            https://hub.docker.com/r/hariszaf/microbetag_prep
        """
        raise ValueError(e)







