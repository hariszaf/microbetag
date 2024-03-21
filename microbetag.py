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
import subprocess
from utils import *
from config import Config
from build_cx_annotated_graph import *
from julia.api import Julia

config_file = sys.argv[1]

with open(config_file, 'r') as yaml_file:
    config = Config(yaml.safe_load(yaml_file), config_file)

if config.bins_path is None:
    raise ValueError

# ----------------
# Build network if not available
# ----------------
if not os.path.exists(config.network) or os.path.getsize(config.network) == 0:
    print("\n >> NETWORK INFERENCE WITH FLASHWEAVE \n")
    ensure_flashweave_format(conf=config)
    pair_args = set()
    for arg, values in config.flashweave_args.items():
        if values["required"]:
            if isinstance(values["value"], bool):
                pair_args.add( ( arg, str(values["value"]).lower()) )
            else:
                print(f'You need to provide values for "{arg}" argument of FlashWeave.') ; sys.exit(0)
        else:
            if values["value"] is not None:
                if values["type"] == "Bool":
                    pair_args.add( (arg, str(values["value"]).lower()) )
                else:
                    pair_args.add( (arg, values["value"]) )

    pair_args.add(("transposed", "true"))
    learn_in = ",".join(f"{arg[0]}={arg[1]}" for arg in pair_args)

    jl = Julia(compiled_modules=False)
    jl.using("FlashWeave")
    if config.metadata_file:
        jl.eval(f'save_network("{config.network}", "{config.metadata_file}", \
            learn_network("{config.flashweave_abd_table}", {learn_in}))')
    else:
        jl.eval(f'save_network("{config.network}", learn_network("{config.flashweave_abd_table}", {learn_in}))')

    ensure_same_namespace_after_fw(config)

# ----------------
# FAPROTAX
# ----------------
print("\n >> LITERATURE ANNOTATION WITH FAPROTAX  \n")
faprotax_params = [
    "python3", config.faprotax_script,
    "-i", config.abundance_table,
    "-o", config.faprotax_funct_table,
    "-g", config.faprotax_txt,
    "-c", '"' + "#" + '"',
    "-d", '"' + config.taxonomy_column_name + '"',
    "-v",
    "--force",
    "-s", config.faprotax_sub_tables,
]
faprotax_command = " ".join(faprotax_params)
process = subprocess.Popen(faprotax_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
stdout, stderr = process.communicate()

# ----------------
# phen annotations
# ----------------
print("\n >> PREDICTING PHENOTYPIC TRAITS \n")
suffixes = [".fa", ".fasta", ".gz"]
bin_files = get_files_with_suffixes(config.bins_path, suffixes)
bin_files_in_a_row = " ".join(bin_files)

# Build genotypes
if get_library_version("scikit-learn") != "1.3.2":
    os.system("python3 -m pip install scikit-learn==1.3.2")
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
    if os.system(compute_genotype_command) != 0:
        print("Try phenotrex genotype for the second time.")
        if os.system(compute_genotype_command) != 0:
            raise ValueError

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
# Prodigal - using DiTing interface
# ----------------
print("\n >> PREDICTING ORFs \n")
# [TODO] Avoid double - prodigal run if user can provide it
for bin_fa in bin_files:
    bin_filename = os.path.basename(bin_fa)
    bin_id, extension = os.path.splitext(bin_filename)
    run_prodigal(bin_fa, bin_id, config.prodigal)

# ----------------
# KEGG annotation - using DiTing interface
# ----------------
print("\n >> KEGG ANNOTATION OF THE PRODIGAL ORFs \n")
ko_list = os.path.join(config.kegg_db_dir, 'ko_list')
ko_dic = ko_list_parser(ko_list)

ko_merged_tab = os.path.join(config.kegg_annotations, 'ko_merged.txt')
if not os.path.exists(ko_merged_tab):
    for bn in config.bin_filenames:
        bin_id, extension = os.path.splitext(bn)
        faa = os.path.join(config.prodigal, bin_id + '.faa')
        kegg_annotation(faa, bin_id, config.kegg_pieces_dir, config.kegg_db_dir, ko_dic, config.threads)

bins_kos, pivot_df = merge_ko(config.kegg_pieces_dir, ko_merged_tab)

# ----------------
# Extract pathway complementarities
# ----------------
print("\n>> GET PATHWAY COMPLEMENTS ")
bin_kos_per_module, alt_to_gapfill, complements = export_pathway_complementarities(
    config,
    pivot_df
)

if not os.path.exists(config.alts_file):
    with open(config.alts_file, "w") as file:
        json.dump(alt_to_gapfill, file, cls=SetEncoder)

if not os.path.exists(config.compl_file):
    with open(config.compl_file, "w") as file:
        json.dump(complements, file, cls=SetEncoder)

# ----------------
# Build GENREs
# ----------------
if config.users_models is False and config.seed_complementarity:
    print("\n >> GENOME-SCALE METABOLIC NETWORK RECONSTRUCTIONS \n")

    build_genres = build_genres(config)

    # Annotate step
    if config.input_for_recon_type == "bins_fasta":

        if config.genre_reconstruction_with == "modelseedpy":
            build_genres.rast_annotate_genomes()
        elif config.gene_predictor == "prodigal":
            print("DiTing .faa files will be used")  # go to the .faa case
        elif config.gene_predictor == "fragGeneScan":
            print("Get annotations with FragGeneScan.")
            build_genres.fgs_annotate_genomes()

    elif config.input_for_recon_type == "coding_regions":
        print("CarveMe will be used with the users .ffn-like files.")

    else:
        print(f"The combination of gene_predictor: {config.gene_predictor} \
            \nand genre_reconstruction_with: {config.genre_reconstruction_with}, are not supported")

    # Reconstruct step
    if config.genre_reconstruction_with == "modelseedpy":
        build_genres.modelseed_reconstructions()
    elif config.genre_reconstruction_with == "carveme":
        build_genres.carve_reconstructions()
    else:
        print("User models to be used for the seed complementarity step.")


# ----------------
# Phylomint
# ----------------
if config.seed_complementarity:
    print("\n >> COMPUTING SEED SETS AND SCORES \n")
    if not os.path.exists(config.phylomint_scores):
        run_phylomint(config)
    else:
        print("Seed scores already computed.")

# ----------------
# Export seed complementarities
# ----------------
if config.seed_complementarity:
    print("\n >> EXPORTING SEED COMPLEMENTS \n")
    seed_complements = export_seed_complementarities(config)
    """
    [NOTE]:consider running again "seed scores" (PhyloMint) using update seed sets
    in this case, we should also edit the ConfidenceScore dictionary
    by removing seeds that were removed in the update()
    """
    if not os.path.exists(seed_complements.updated_seed_sets):
        print("!!! Updating seed and non seed sets !!! ")
        seed_complements.update()
    else:
        print("Seed sets already updated.")

    if config.genre_reconstruction_with == "carveme":
        print("We will map the BIGG compounds to ModelSEED ones.\
            \nIn the future, we will map BiGG ids to KEGG so we do not have to go through ModelSEED in this scenario.")
        seed_complements.map_carveme_seeds()

    if not os.path.exists(seed_complements.module_seeds):
        seed_complements.module_related_seeds()
    else:
        print("Seed and non seed sets with compounds related to KEGG modules already retrieved.")

    if not os.path.exists(seed_complements.seed_complements):
        seed_complements.export_seed_complements()
    else:
        print("Seed complements already exported.")


# ----------------
# Annotate network in .cx format
# ----------------
print("\n >> ANNOTATE NETWORK \n")
annotated_network = build_cx_annotated_graph(config)
with open(config.microbetag_annotated_network_file, "w") as f:
        annotated_network2file = convert_to_json_serializable(annotated_network)
        json.dump(annotated_network2file, f)
