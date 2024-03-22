
import os
import pandas as pd
import sys
import cobra

class Config:
    """
    Parses a yaml file to init values.
    """
    def __init__(self, conf, config_file):

        # User's group and user id
        file_info = os.stat(config_file)
        self.user_id = file_info.st_uid
        self.group_id = file_info.st_gid

        self.io_path = conf["io_path"]["path"]
        self.cwd = os.getcwd()
        self.kegg_db_dir = os.path.join(self.cwd, "microbetagDB/ref-dbs/kofam_database/")
        self.threads = conf["threads"]["value"]

        self.on_container = True if conf["on_container"]["value"] else False
        self.mount = "/data" if conf["on_container"]["value"] else self.io_path
        self.bins_path = os.path.join(self.mount, conf["bins_fasta"]["folderName"])
        self.abundance_table = os.path.join(self.mount, conf["abundance_table_file"]["fileName"])
        self.output_dir = os.path.join(self.mount, conf["output_directory"]["folderName"])
        self.network = os.path.join(self.mount, conf["edge_list"]["value"]) if conf["edge_list"]["value"] else None
        self.metadata_file = os.path.join(self.mount, conf["metadata_file"]["fileName"]) if conf["metadata_file"]["fileName"] and os.path.exists(os.path.join(self.mount, conf["metadata_file"]["fileName"])) else None
        self.flashweave_abd_table = os.path.join(self.mount, "abd_table_for_flashweave.tsv")

        self.bin_filenames = os.listdir(self.bins_path)

        self.input_for_recon_type = conf["input_type_for_seed_complementarities"]["value"] \
            if conf["input_type_for_seed_complementarities"]["value"] in conf["input_type_for_seed_complementarities"]["value_from"] \
            else sys.exit("Invalid input_type_for_seed_complementarities specified.")

        self.seed_complementarity = conf["seed_complementarity"]["value"]

        input_value = conf["input_type_for_seed_complementarities"]["value"]
        allowed_values = conf["input_type_for_seed_complementarities"]["value_from"]
        self.input_for_recon_type = input_value if input_value in allowed_values else (print("Error: Input value is not among the allowed values:", allowed_values) or sys.exit(1))

        self.users_models = True if conf["input_type_for_seed_complementarities"]["value"] == "models" else False
        self.for_reconstructions = os.path.join(self.mount, conf["sequence_files_for_reconstructions"]["folderName"])

        # Check whethere bin names are the same in both abundance and edgelist files
        bins = [ os.path.splitext(gbin)[0] for gbin in self.bin_filenames  ]
        f = pd.read_csv(self.abundance_table, sep="\t")
        self.taxonomy_column_name = list(f.columns)[-1]
        self.sequence_column_name = list(f.columns)[0]
        bins_in_abundance_file = list(f.iloc[:,0])

        if not all(elem in bins_in_abundance_file for elem in bins):
            not_in_second_list = set(bins) - set(bins_in_abundance_file)
            not_in_second_list_str = ', '.join(not_in_second_list)
            raise ValueError(f"Bin names do not match with those in the abundance table: {not_in_second_list_str}")

        if self.network:
            self.flashweave = False
            f = pd.read_csv(self.network, sep="\t")
            bins_in_net = set(list(f.iloc[:,0]) + list(f.iloc[:,1]))
            if not all(elem in bins_in_abundance_file for elem in bins_in_net) or not all(elem in bins_in_net for elem in bins):
                raise ValueError(f"Bin names in the edgelist file do not match with those in the abundance table and/or in the bins.")
        else:
            self.flashweave_script = os.path.join(self.cwd, "microbetagDB/scripts/flashweave.jl")
            self.network = os.path.join(self.output_dir, "network_output.edgelist")
            self.flashweave = True

        # Build output dir
        os.makedirs(self.output_dir, exist_ok=True)

        self.predictions_path = os.path.join(self.output_dir, "predictions")
        os.makedirs(self.predictions_path, exist_ok=True)

        self.prodigal = os.path.join(self.output_dir, "ORFs")
        os.makedirs(self.prodigal, exist_ok=True)

        self.kegg_annotations = os.path.join(self.output_dir, "KEGG_annotations")
        self.kegg_pieces_dir = os.path.join(self.kegg_annotations, 'hmmout')
        os.makedirs(self.kegg_annotations, exist_ok=True)
        os.makedirs(self.kegg_pieces_dir, exist_ok=True)


        self.reconstructions = os.path.join(self.output_dir, "reconstructions")
        self.genres = os.path.join(self.reconstructions, "GENREs")
        os.makedirs(self.reconstructions, exist_ok=True)
        os.makedirs(self.genres, exist_ok=True)

        self.gene_predictor = conf["gene_predictor"]["value"]
        self.genre_reconstruction_with = conf["genre_reconstruction_with"]["value"]

        self.seeds = os.path.join(self.output_dir, "seeds_complementarity")
        os.makedirs(self.seeds, exist_ok=True)
        self.seed_complements = os.path.join(self.seeds, "seed_complements.pckl")
        self.module_related_non_seeds = os.path.join(self.seeds, "module_related_non_seeds.pckl")
        self.phylomint_scores = os.path.join(self.seeds, "phylomint_scores.tsv")

        self.pathway_complements_dir = os.path.join(self.output_dir, "pathway_complementarity")
        os.makedirs(self.pathway_complements_dir, exist_ok=True)
        self.alts_file = os.path.join(self.pathway_complements_dir, "alts.json")
        self.compl_file = os.path.join(self.pathway_complements_dir, "pathCompls.json")
        self.pathway_complement_percentage = conf["pathway_complement_percentage"]["value"] if conf["pathway_complement_percentage"]["value"] is not None else 0

        # ModelSEEDpy arguments
        self.gapfill_model = conf["gapfill_model"]["value"]
        self.gapfill_media = conf["gapfill_media"]["value"]

        # Flashweave arguments
        self.metadata = "false" if self.metadata_file == "false" else "true"
        self.flashweave_args = conf["flashweave_args"]

        # Phenotrex
        self.genotypes_file = os.path.join(self.output_dir, "train.genotype")
        self.min_proba = conf["min_proba"]["value"]

        # Mappings
        self.kegg_mappings = os.path.join(self.cwd, "microbetagDB/mappings/kegg_mappings/")
        self.ko_terms_per_module_definition = os.path.join(self.kegg_mappings, "kegg_terms_per_module.tsv")
        self.modules_definitions_json_map = os.path.join(self.kegg_mappings, "module_definition_map.json")
        self.kegg_modules_to_maps = os.path.join(self.kegg_mappings, "module_map_pairs.tsv")
        self.seed_ko_mo = os.path.join(self.kegg_mappings, "seedId_keggId_module.tsv")
        self.module_descriptions = os.path.join(self.kegg_mappings, "module_descriptions")
        self.metanetx_compounds = os.path.join(self.cwd, "microbetagDB/mappings/MetaNetX/chem_xref.tsv")

        # FAPROTAX
        self.faprotax_txt = os.path.join(self.cwd, "microbetagDB/ref-dbs/FAPROTAX_1.2.7/FAPROTAX.txt")
        self.faprotax_script = os.path.join(self.cwd, "microbetagDB/ref-dbs/FAPROTAX_1.2.7/collapse_table.py")
        self.faprotax_output_dir = os.path.join(self.output_dir, "faprotax")
        self.faprotax_funct_table = os.path.join(self.faprotax_output_dir, "functional_otu_table.tsv")
        self.faprotax_sub_tables = os.path.join(self.faprotax_output_dir, "sub_tables")
        os.makedirs(self.faprotax_output_dir, exist_ok=True)
        os.makedirs(self.faprotax_sub_tables, exist_ok=True)


        self.network_clustering = conf["network_clustering"]["value"] if conf["network_clustering"]["value"] else False

        self.max_scratch_alt = conf["max_length_for_complement_from_scratch"]["value"] if conf["max_length_for_complement_from_scratch"]["value"] else 1

        self.microbetag_annotated_network_file = os.path.join(self.output_dir, "microbetag_annotated_network.cx")


        # ==========
        # Init torch
        # ==========
        import torch
        from deepnog.utils import get_weights_path
        from deepnog.utils import set_device
        device = set_device('auto')
        weights_path = get_weights_path(
            database="eggNOG5",
            level=str(2),
            architecture="deepencoding",
        )
        model_dict = torch.load(weights_path, map_location=device)

        # Tests
        if self.users_models:
            random_model = os.path.join(self.for_reconstructions, os.listdir(self.for_reconstructions)[0])
            model = cobra.io.read_sbml_model(random_model)
            if model.metabolites[0].id[:3] == "cpd" and self.genre_reconstruction_with == "carveme":
                print("Based on the genre_reconstruction_with argument, your model are expected to have been built using the BiGG namespace\
                    \nyet, they are using the ModelSEED one. Please make sure you set those arguments in line or use another set of models that use the BiGG namespace indeed.");sys.exit(0)
            elif model.metabolites[0].id[:3] == "cpd" and self.genre_reconstruction_with != "modelseedpy":
                self.genre_reconstruction_with = "modelseedpy"
                print("WARNING. Var was reset")
            elif model.metabolites[0].id[:3] != "cpd" and self.genre_reconstruction_with == "modelseedpy":
                print("Based on the genre_reconstruction_with argument, your model are expected to have been built\
                      \nusing ModelSEED but their namespace does not aggre. Make sure the genre_reconstruction_with variable aggree with your models format.");sys.exit(0)
            elif model.metabolites[0].id[:3] != "cpd" and self.genre_reconstruction_with != "carveme":
                print("WARNING! The models are assumed to use the BiGG namespace.")
                self.genre_reconstruction_with = "carveme"