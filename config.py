
import os

class Config:
    """
    Parses a yaml file to init values.
    """
    def __init__(self, conf):

        # # User's group and user id
        # file_info = os.stat(config_file)
        # self.user_id = file_info.st_uid
        # self.group_id = file_info.st_gid

        self.io_path = conf["io_path"]["path"]

        if conf["on_container"]["value"]:
            self.mount = "/data"
            self.bins_path = os.path.join(self.mount, conf["bins_fasta"]["folderName"])
            self.abundance_table = os.path.join(self.mount, conf["abundance_table_file"]["fileName"])
            self.output_dir = os.path.join(self.mount, conf["output_directory"]["folderName"])
        else:
            self.bins_path = os.path.join(self.io_path, conf["bins_fasta"]["folderName"])
            self.abundance_table = os.path.join(self.io_path, conf["abundance_table_file"]["fileName"])
            self.output_dir = os.path.join(self.io_path, conf["output_directory"]["folderName"])


        if os.path.exists(self.output_dir) and os.path.isdir(self.output_dir):
            print("Output directory already exists.")
        else:
            os.mkdir(self.output_dir)

        self.predictions_path = os.path.join(self.output_dir, "predictions")
        if os.path.exists(self.predictions_path) and os.path.isdir(self.predictions_path):
            print("Predictions path already exists.")
        else:
            os.mkdir(self.predictions_path)

        # Metadata file for FlashWeave
        if conf["metadata_file"]["fileName"]:
            self.metadata = "true"
            if conf["on_container"]["value"]:
                self.metadata_file = os.path.join(self.mount, conf["metadata_file"]["fileName"])
            else:
                self.metadata_file = os.path.join(self.io_path, conf["metadata_file"]["fileName"])
        else:
            self.metadata = "false"


        self.build_network = conf["build_network"]["value"]
        self.threads = conf["threads"]["value"]

        # Flashweave arguments
        if conf["flashweave_sensitive"]["value"]:
            self.sensitive = "true"
        else:
            self.sensitive = "false"
        if conf["flashweave_heterogeneous"]["value"]:
            self.heterogeneous = "true"
        else:
            self.heterogeneous = "false"


        self.genotypes_file = os.path.join(self.output_dir, "train.genotype")
        self.min_proba = conf["min_proba"]["value"]



        # init torch
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

