# microbetag : a software suite to annotate microbial co-occurrence networks

# Copyright (c) 2025 Haris Zafeiropoulos

# Licensed under GNU LGPL.3, see LICENCE file

import os
import json
import pandas as pd

from .helpers import (
    Emojis,
    Faprotax,
    MappingPaths,
    NetworkHandler,
    AbdTableHandler,
    PathwayComplementarity,
    SeedComplementarityHandler,
    missing_input_error,
)
from .networks import get_edgelist
from .utils import (
    mtg_logger,
    detect_separator,
    resolve_file_path,
    convert_to_json_serializable
)


_logger_ = mtg_logger(__name__)


class Config:
    """
    Parses a microbetag configuration file (yaml) to init a microbetag run.

    Args:
        conf: A dictionary where the YAML configuration file has been loaded
        config_file: Filepath to the configuration YAML file.

    Attention:
        Use the configuration template of the microbetag version you are running
        (`microbetag --template config.yml`); otherwise the Config class may fail to create an instance.

    Example:
        >>> with open(args.config, "r") as yaml_file:
                yaml_conf = yaml.safe_load(yaml_file)
        >>> conf = Config(yaml_conf, args.config)
    """

    def __init__(self, conf: dict, config_file: str = None):

        # Pass loaded yaml object
        self.yaml = conf

        # Load emojis
        emojis = Emojis()

        # User's group and user id
        if config_file is not None:
            self.yaml_file = os.stat(config_file)
            self.user_id   = self.yaml_file.st_uid
            self.group_id  = self.yaml_file.st_gid

        self.cwd = os.path.dirname(os.path.realpath(__file__))

        if config_file is None:
            self.base_dir = conf.get("__config_dir__", None)
        else:
            self.base_dir = os.path.dirname(config_file)
        if self.base_dir is None:
            raise SystemError("You need to provide a base directory, where your config is found.")

        # Output dir
        outdir = conf.get("output_directory", {}).get("dir_path")
        if outdir:
            self.outdir = os.path.join(self.base_dir, outdir)
        else:
            raise ValueError("Output directory needs to be specified.")

        # Build output dir
        os.makedirs(self.outdir, exist_ok=True)

        # Mappings
        mappings = MappingPaths()
        self.__dict__.update(vars(mappings))

        # Checks whether microbetag is running on the on-the-fly version, by default False
        self.onthefly = get_value(conf, "onthefly", False)
        self.api      = get_value(conf, "api", False)

        # Threads to be used
        self.threads = get_value(conf, "threads", 2)

        # Steps
        self.faprotax    = get_value(conf, "faprotax_annotation", False)
        self.phen_traits = get_value(conf, "phenotrex_traits", False)
        self.path_compl  = get_value(conf, "pathway_complementarity", False)
        self.seed_compl  = get_value(conf, "seed_complementarity", False)
        self.net_cluster = get_value(conf, "network_clustering", False)

        # The abundance table is optional; without it, a sequence id to taxonomy map is needed
        abd_tbl_filename     = conf.get("abundance_table_file", {}).get("file_path")
        self.abundance_table = resolve_file_path(self.base_dir, abd_tbl_filename)

        # Edgelist of provided network
        edge_list    = conf.get("edge_list", {}).get("file_path")
        self.network = resolve_file_path(
            self.base_dir, edge_list
        )

        # -------
        # NOTE: Setting precalculations only as true allows to get a Config instance
        # without providing an abundance table or network, meaning you can use the Config instance
        # for partial/specific tasks of microbetag.
        # -------

        precalc_only      = get_value(conf, "precalculations_only", False)
        self.precalc_only = precalc_only if precalc_only in [0, 1] else False

        if self.network is None and not self.precalc_only:
            raise ValueError(
                "Please provide the co-occurrence network to annotate (`edge_list: file_path`), "
                "as a 3-column edge list. To infer one from your abundance table with FlashWeave, "
                "use the microbetag-nf pipeline (https://github.com/msysbio/microbetag-nf)."
            )

        # IMPORTANT: Sequence to taxonomy map -- required if no abundance table is needed
        sequence_taxonomy_map      = conf.get("sequence_id_taxonomy_map", {}).get("file_path")
        self.sequence_taxonomy_map = resolve_file_path(
            self.base_dir, sequence_taxonomy_map
        )

        if (
            self.abundance_table is None and self.sequence_taxonomy_map is None and not self.precalc_only
        ):
            raise ValueError(
                "Since an abundance table is not provided, you need to provide a 2-column file "
                "with the sequence id (e.g bin ids) and their corresponding taxonomy or taxon name."
                "Thus, you need to include the `sequence_id_taxonomy_map` parameter on your config file."
                "Check the complete configuration file for that."
            )

        # IMPORTANT: Sequence id to taxonomy map
        if self.network is None and self.abundance_table:
            (
                self.seq_to_taxon_df,
                self.sequence_id_column_name,
                self.taxonomy_column_name,
                self.delimiter,

            ) = load_abundance(self.abundance_table)

            self.seq_ids = self.seq_to_taxon_df["sequence_id"].unique().tolist()

        elif self.abundance_table is None and self.network:

            self.delimiter = detect_separator(self.sequence_taxonomy_map)

            seq_to_taxon_df         = pd.read_csv(self.sequence_taxonomy_map, sep=self.delimiter)
            seq_to_taxon_df.columns = ["sequence_id", "taxonomy"]
            self.seq_to_taxon_df    = seq_to_taxon_df
            self.seq_ids            = self.seq_to_taxon_df["sequence_id"].unique().tolist()

        elif self.abundance_table and self.network:

            # -------
            # NOTE: Not all sequence ids in the seq_ids need to have a taxonomy in this case --
            # only those coming from the abundance table
            # Yet, in case that the network has taxa not present in the abundance table, it will lead to errors.
            # -------

            network_df  = get_edgelist(self.network)
            net_seq_ids = (
                pd.concat([network_df.iloc[:, 0], network_df.iloc[:, 1]])
                .unique()
                .tolist()
            )

            (
                self.seq_to_taxon_df,
                self.sequence_id_column_name,
                self.taxonomy_column_name,
                self.delimiter,

            ) = load_abundance(self.abundance_table)

            abd_seq_ids = self.seq_to_taxon_df["sequence_id"].unique().tolist()
            self.seq_ids = net_seq_ids + abd_seq_ids

        else:
            _logger_.warning(
                f"{emojis.WARNING_EMOJI}Neither an abundance table nor a network was procided.\n"
                "microbetag will only run some pre-calculations not requiring them.\n"
            )

        # Load abundance table with taxonomy
        nsc = AbdTableHandler(config=self)
        self.__dict__.update(vars(nsc))

        # Check whether bin names are the same in both abundance and edgelist files
        nh = NetworkHandler(config=self)
        self.__dict__.update(vars(nh))

        # Update conf variables for pathway complementarity module
        if self.path_compl:
            pc = PathwayComplementarity(config=self)
            self.__dict__.update(vars(pc))

        # Phenotrex: the predictions are an input; in the on-the-fly version they come from microbetagDB
        if self.phen_traits:
            if self.onthefly:
                self.phen_traits_dir = os.path.join(self.outdir, "phen_predictions")
                os.makedirs(self.phen_traits_dir, exist_ok=True)
            else:
                predictions = conf.get("phenotrex_predictions", {}).get("dir_path")
                if predictions is None:
                    raise missing_input_error(
                        "Phenotypic trait annotation", "phenotrex_predictions", "folder of predictions",
                        "phenotrex (phenotrex predict, one .tsv per trait)",
                    )
                self.phen_traits_dir = resolve_file_path(self.base_dir, predictions)

        # FAPROTAX
        if self.abundance_table is not None and self.faprotax:
            faprotax = Faprotax(config=self)
            self.__dict__.update(vars(faprotax))

        # manta: the clustered network is an input
        if self.net_cluster:
            manta_net = conf.get("prev_clustered_network", {}).get("file_path")
            if manta_net is None:
                raise missing_input_error(
                    "Network clustering", "prev_clustered_network", "clustered network (.cyjs)",
                    "manta (its input can be built with microbetag-format-manta)",
                )
            self.manta_net = resolve_file_path(self.base_dir, manta_net)

        # Seed complementarity
        if self.seed_compl:
            sc = SeedComplementarityHandler(config=self)
            self.__dict__.update(vars(sc))

        # Intermediate annoteted network file name
        self.microbetag_annotated_network_file = os.path.join(
            self.outdir, "pseudo_cx_annotated_net.cx"
        )
        self.tinyurl = get_value(conf, "tinyurl", False)

        _logger_.info("Configuration file loaded successfully.")

    def export_to_log(self, log_file="parameters.log"):
        """Dumps the Config instance in a JSON file."""
        args = convert_to_json_serializable(self.__dict__)
        with open(log_file, "w") as f:
            json.dump(args, f)


def get_value(conf, key, default=None):
    """
    Retrieves the 'value' field for a given key from a nested configuration dictionary.

    Returns:
    any: The value associated with conf[key]['value'], or the provided default if not found or None.
    """
    value = conf.get(key, {}).get("value", default)
    return default if value is None else value


def load_config(yaml_file):

    import yaml

    config_path = os.path.abspath(yaml_file)
    config_dir  = os.path.dirname(config_path)

    with open(yaml_file, "r") as y:
        yaml_data = yaml.safe_load(y)
    yaml_data["__config_dir__"] = config_dir

    return yaml_data


def load_abundance(abd_file: str) -> tuple[pd.DataFrame, str, str, str]:
    """
    Load a tsv/csv format abundance table assuming the sequence id is procided in the first column
    and the taxonomy in the last one

    Args:
        abd_file: Filepath to abundance table file.

    Returns:
        A tuple including:
            - seq_id2tax: A :class:`pandas.DataFrame` with the sequence id and their corresponding taxonomy
            - seq_id_col: The name of the column with the sequence identifier (e.g. ``seqId``)
            - tax_col: The name of the column with the taxonomy
    """

    delimiter          = detect_separator(abd_file)
    abd_tab_df         = pd.read_csv(abd_file, sep=delimiter)
    seq_id_col         = abd_tab_df.columns[0]
    tax_col            = abd_tab_df.columns[-1]
    seq_id2tax         = abd_tab_df[[seq_id_col, tax_col]]
    seq_id2tax.columns = ["sequence_id", "taxonomy"]

    return seq_id2tax, seq_id_col, tax_col, delimiter
