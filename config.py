
import os
import pandas as pd

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

        self.mount = "/data" if conf["on_container"]["value"] else self.io_path
        self.bins_path = os.path.join(self.mount, conf["bins_fasta"]["folderName"])
        self.abundance_table = os.path.join(self.mount, conf["abundance_table_file"]["fileName"])
        self.output_dir = os.path.join(self.mount, conf["output_directory"]["folderName"])
        self.network = os.path.join(self.mount, conf["edge_list"]["value"]) if conf["edge_list"]["value"] else None
        self.metadata_file = os.path.join(self.mount, conf["metadata_file"]["fileName"]) if conf["metadata_file"]["fileName"] else "false"
        self.flashweave_abd_table = os.path.join(self.mount, "abd_table_for_flashweave.tsv")

        self.bin_filenames = os.listdir(self.bins_path)


        # Check whethere bin names are the same in both abundance and edgelist files
        bins = [ os.path.splitext(gbin)[0] for gbin in self.bin_filenames  ]
        f = pd.read_csv(self.abundance_table, sep="\t")
        bins_in_abundance_file = list(f.iloc[:,0])

        if not all(elem in bins_in_abundance_file for elem in bins):
            not_in_second_list = set(bins) - set(bins_in_abundance_file)
            not_in_second_list_str = ', '.join(not_in_second_list)
            raise ValueError(f"Bin names do not match with those in the abundance table: {not_in_second_list_str}")

        if self.network:
            f = pd.read_csv(self.network, sep="\t")
            bins_in_net = set(list(f.iloc[:,0]) + list(f.iloc[:,1]))
            if not all(elem in bins_in_abundance_file for elem in bins_in_net) or not all(elem in bins_in_net for elem in bins):
                raise ValueError(f"Bin names in the edgelist file do not match with those in the abundance table and/or in the bins.")
        else:
            self.flashweave_script = os.path.join(self.cwd, "microbetagDB/scripts/flashweave.jl")


        # Build output dir
        os.makedirs(self.output_dir, exist_ok=True)

        self.predictions_path = os.path.join(self.output_dir, "predictions")
        os.makedirs(self.predictions_path, exist_ok=True)

        self.prodigal = os.path.join(self.output_dir, "ORFs")
        os.makedirs(self.prodigal, exist_ok=True)

        self.reconstructions = os.path.join(self.output_dir, "reconstructions")
        self.genres = os.path.join(self.reconstructions, "GENREs")
        os.makedirs(self.reconstructions, exist_ok=True)
        os.makedirs(self.genres, exist_ok=True)

        self.kegg_annotations = os.path.join(self.output_dir, "KEGG_annotations")
        self.kegg_pieces_dir = os.path.join(self.kegg_annotations, 'hmmout')
        os.makedirs(self.kegg_annotations, exist_ok=True)
        os.makedirs(self.kegg_pieces_dir, exist_ok=True)

        self.seeds = os.path.join(self.output_dir, "seeds")
        os.makedirs(self.seeds, exist_ok=True)

        # ModelSEEDpy arguments
        self.gapfill_model = conf["gapfill_model"]["value"]
        self.gapfill_media = conf["gapfill_media"]["value"]

        # Flashweave arguments
        self.build_network = conf["build_network"]["value"]
        self.metadata = "false" if self.metadata_file == "false" else "true"
        self.sensitive = "true" if conf["flashweave_sensitive"]["value"] else "false"
        self.heterogeneous = "true" if conf["flashweave_heterogeneous"]["value"] else "false"

        # Phenotrex
        self.genotypes_file = os.path.join(self.output_dir, "train.genotype")
        self.min_proba = conf["min_proba"]["value"]

        self.phylomint_scores = os.path.join(self.seeds, "phylomint_scores.tsv")

        self.ko_terms_per_module_definition = os.path.join(self.cwd, "microbetagDB/mappings/kegg_mappings/kegg_terms_per_module.tsv")
        self.modules_definitions_json_map = os.path.join(self.cwd, "microbetagDB/mappings/kegg_mappings/module_definition_map.json")
        self.kegg_modules_to_maps = os.path.join(self.cwd, "microbetagDB/mappings/kegg_mappings/module_map_pairs.tsv")
        self.seed_ko_mo = os.path.join(self.cwd, "microbetagDB/mappings/kegg_mappings/seedId_keggId_module.tsv")

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

