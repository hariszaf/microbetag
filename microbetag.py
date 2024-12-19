"""
Aim:
    Perform the microbetag () approach annotating user's bins instead of mapping taxonomies against the
    microbetagDB representative genomes.

Input:
    - A folder with .fasta files of the corresponding bins
    - The abundance table of the bins across the samples
    - (optional) a co-occurrence network in a three-columns format

Output:
    - An annotated network in .cx format

Author:
    Haris Zafeiropoulos

"""

__version__ = "v1.0.2"

import os
import sys
import logging

# Set up custom logging format
logging.basicConfig(
    format='%(levelname)s: %(message)s',  # Define the format without "root:"
    level=logging.INFO  # Set the logging level
)


def print_help():
    help_message = """
    Usage: python microbetag.py <path_to_config_yml>

    Other options:
      h        Display this help message.
      v        Display version.
    """
    print(help_message)


def print_version():
    print(__version__)

def print_config_message():
    conf_message = """
    The config file you provided cannot be loaded.
    Please make sure you follow the instructions on the documentation site:
    https://hariszaf.github.io/microbetag/docs/tutorials/local/#input-and-configyml-files
    """
    logging.error(conf_message)


if len(sys.argv) == 1:
    print_help(); sys.exit()

if sys.argv[1] == '-h' or sys.argv[1] == '--help':
    print_help(); sys.exit()

if sys.argv[1] == 'v' or sys.argv[1] == 'version':
    print_version(); sys.exit()

import yaml
import subprocess
from utils import *
from config import Config
from build_cx_annotated_graph import *
from julia.api import Julia

config_file = sys.argv[1]

with open(config_file, 'r') as yaml_file:
    try:
        config = Config(yaml.safe_load(yaml_file), config_file)
    except:
        print_config_message() ; sys.exit(0)

if config.bins_path is None:
    raise ValueError

# ----------------
# Build network if not available
# ----------------
if config.precalc_only:
    logging.warning("No abundance table or network was provided. microbetag will try to run pre-calculations.")
elif not os.path.exists(config.network) or os.path.getsize(config.network) == 0:
    logging.info("\n >> NETWORK INFERENCE WITH FLASHWEAVE \n")
    ensure_flashweave_format(conf=config)
    pair_args = set()
    for arg, values in config.flashweave_args.items():
        if values["required"]:
            if isinstance(values["value"], bool):
                pair_args.add( ( arg, str(values["value"]).lower()) )
            else:
                logging.error(f'You need to provide values for "{arg}" argument of FlashWeave.') ; sys.exit(0)
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
if config.abundance_table is not None:
    logging.info("[STEP] LITERATURE ANNOTATION WITH FAPROTAX")
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
    logging.info("[STEP] PREDICTING PHENOTYPIC TRAITS")
    if os.system(compute_genotype_command) != 0:
        logging.info("Try phenotrex genotype for the second time.")
        if os.system(compute_genotype_command) != 0:
            logging.error("asda")
            sys.exit(0)

# Get predictions
folder_path = "microbetagDB/ref-dbs/phenDB/classes/"
phen_models = [os.path.join(folder_path, model) for model in os.listdir(folder_path)]

for model in phen_models:
    model_name =  os.path.basename(model)
    model_predictions_output = "".join([
        config.predictions_path, "/", model_name[:-4], ".prediction.tsv"
    ])
    if os.path.exists(model_predictions_output):
        logging.info("Predictions already exist for model: %s", model_name)
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
            logging.error("TSIRIMPIM") ; sys.exit(0)

# ----------------
# Prodigal - using DiTing interface
# ----------------
logging.info("[STEP  ] PREDICTING ORFs WITH PRODIGAL THROUGH DiTing")
# [TODO] Avoid double - prodigal run if user can provide it
for bin_fa in bin_files:
    bin_filename = os.path.basename(bin_fa)
    bin_id, extension = os.path.splitext(bin_filename)
    run_prodigal(bin_fa, bin_id, config.prodigal)



# ----------------
# Pathway complementarity
# ----------------
if config.pathway_complementarity:
    # ----------------
    # KEGG annotation - using DiTing interface // required in case of pathway complementarities
    # ----------------
    ko_list = os.path.join(config.kegg_db_dir, 'ko_list')
    ko_dic = ko_list_parser(ko_list)

    if config.ko_merged is None:
        config.ko_merged = os.path.join(config.kegg_annotations, 'ko_merged.txt')
        logging.info("[STEP ] KEGG ANNOTATION OF THE PRODIGAL ORFs \n")
        for bn in config.bin_filenames:
            bin_id, extension = os.path.splitext(bn)
            faa = os.path.join(config.prodigal, bin_id + '.faa')
            kegg_annotation(faa, bin_id, config.kegg_pieces_dir, config.kegg_db_dir, ko_dic, config.threads)

        merge_ko(config.kegg_pieces_dir, config.ko_merged)
    else:
        logging.info("A 3-col KEGG annotation file already available.")

    pivot_df = load_merged_ko_file(config.ko_merged)  # Load ko_merged.txt

    # ----------------
    # Extract pathway complementarities
    # ----------------
    if not os.path.exists(config.alts_file) or not os.path.exists(config.compl_file):

        logging.info("[STEP ] GET PATHWAY COMPLEMENTS ")
        bin_kos_per_module, alt_to_gapfill, complements = export_pathway_complementarities(
            config,
            pivot_df
        )

# ----------------
# Build GENREs
# ----------------
if config.users_models is False and config.seed_complementarity:

    logging.info("[STEP] GENOME-SCALE METABOLIC NETWORK RECONSTRUCTIONS")

    build_genres = build_genres(config)

    # Annotate step
    if config.input_for_recon_type == "bins_fasta":

        if config.genre_reconstruction_with == "modelseedpy":
            build_genres.rast_annotate_genomes()  # saves under config.reconstructions

        elif config.gene_predictor == "prodigal":
            logging.info("DiTing .faa files will be used")  # go to the .faa case, i.e., the ORFs/

        elif config.gene_predictor == "fragGeneScan":
            logging.info("Get annotations with FragGeneScan.")
            build_genres.fgs_annotate_genomes()   # saves under config.reconstructions

    elif config.input_for_recon_type == "coding_regions":
        logging.info("CarveMe will be used with the users .ffn-like files.")

    else:
        logging.warning(f"The combination of gene_predictor: {config.gene_predictor} \
            \nand genre_reconstruction_with: {config.genre_reconstruction_with}, are not supported")

    # Reconstruct step
    if config.genre_reconstruction_with == "modelseedpy":
        build_genres.modelseed_reconstructions()

    elif config.genre_reconstruction_with == "carveme":
        build_genres.carve_reconstructions()

    else:
        logging.info("User models to be used for the seed complementarity step.")


# ----------------
# Phylomint
# ----------------
if config.seed_complementarity:
    logging.info("[STEP] COMPUTING SEED SETS AND SCORES")
    if not os.path.exists(config.phylomint_scores):
        run_phylomint(config)
    else:
        logging.info("Seed scores already computed.")

# ----------------
# Export seed complementarities
# ----------------
if config.seed_complementarity:
    logging.info("[STEP] EXPORTING SEED COMPLEMENTS")
    seed_complements = export_seed_complementarities(config)
    """
    [NOTE]:consider running again "seed scores" (PhyloMint) using update seed sets
    in this case, we should also edit the ConfidenceScore dictionary
    by removing seeds that were removed in the update()
    """
    if not os.path.exists(seed_complements.updated_seed_sets):
        logging.info("Updating seed and non seed sets!")
        seed_complements.update()
    else:
        logging.info("Seed sets already updated.")

    if config.genre_reconstruction_with == "carveme":
        logging.info("We will map the BIGG compounds to ModelSEED ones.\
            \nIn the future, we will map BiGG ids to KEGG so we do not have to go through ModelSEED in this scenario.")
        seed_complements.map_carveme_seeds()

    if not os.path.exists(seed_complements.module_seeds):
        seed_complements.module_related_seeds()
    else:
        logging.info("Seed and non seed sets with compounds related to KEGG modules already retrieved.")

    if not os.path.exists(seed_complements.seed_complements):
        seed_complements.export_seed_complements()
        logging.info("Seed complements were exported fine.")
    else:
        logging.info("Seed complements already exported.")


# ----------------
# Annotate network in .cx format
# ----------------
if config.precalc_only is False:
    logging.info("[STEP] ANNOTATE NETWORK ")
    annotated_network = build_cx_annotated_graph(config)
    with open(config.microbetag_annotated_network_file, "w") as f:
            annotated_network2file = convert_to_json_serializable(annotated_network)
            json.dump(annotated_network2file, f)
    logging.info("A microbetag-annotated network in .cx format was built sucessfully.")

logging.info("microbetag completed.")
