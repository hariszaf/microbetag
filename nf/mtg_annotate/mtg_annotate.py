import sys
import yaml
import pandas as pd
from pathlib import Path
from microbetag.config import load_abundance
from microbetag.networks import get_edgelist
from microbetag.utils import detect_separator
from microbetag.build_mtg_cx2 import mtg_annotate_network


class Config:
    def __init__(self, yaml_path: str):
        with open(yaml_path, "r") as f:
            data = yaml.safe_load(f)

        if not isinstance(data, dict):
            raise ValueError(f"YAML file {yaml_path} does not contain a valid mapping.")

        for key, value in data.items():
            setattr(self, key, value)

        self.kegg_mappings       = Path("/workspace/microbetag/mtg_maps_models/kegg_mappings/")
        self.module_descriptions = self.kegg_mappings / "module_descriptions"

    def __repr__(self):
        # __repr__ is a special method (a dunder method, short for double underscore) that defines how an object is represented as a string — primarily for developers (not end users).
        # It’s what you see when you type an object’s name in a Python shell or print it inside a list/dict.
        attrs = ", ".join(f"{k}={v!r}" for k, v in self.__dict__.items())
        return f"{self.__class__.__name__}({attrs})"

c = Config(sys.argv[1])

# IMPORTANT: Sequence id to taxonomy map
if c.network is None and c.abundance_file:
    (
        c.seq_to_taxon_df,
        c.sequence_id_column_name,
        c.taxonomy_column_name,
        c.delimiter,

    ) = load_abundance(c.abundance_file)

    c.seq_ids = c.seq_to_taxon_df["sequence_id"].unique().tolist()

elif c.abundance_file is None and c.network:

    c.delimiter = detect_separator(c.sequence_taxonomy_map)

    seq_to_taxon_df         = pd.read_csv(c.sequence_taxonomy_map, sep=c.delimiter)
    seq_to_taxon_df.columns = ["sequence_id", "taxonomy"]
    c.seq_to_taxon_df    = seq_to_taxon_df
    c.seq_ids            = c.seq_to_taxon_df["sequence_id"].unique().tolist()

elif c.abundance_file and c.network:

    # -------
    # NOTE: Not all sequence ids in the seq_ids need to have a taxonomy in this case --
    # only those coming from the abundance table
    # Yet, in case that the network has taxa not present in the abundance table, it will lead to errors.
    # -------

    network_df  = get_edgelist(c.network)
    net_seq_ids = (
        pd.concat([network_df.iloc[:, 0], network_df.iloc[:, 1]])
        .unique()
        .tolist()
    )

    (
        c.seq_to_taxon_df,
        c.sequence_id_column_name,
        c.taxonomy_column_name,
        c.delimiter,

    ) = load_abundance(c.abundance_file)

    abd_seq_ids = c.seq_to_taxon_df["sequence_id"].unique().tolist()
    c.seq_ids   = net_seq_ids + abd_seq_ids

# Annotate
mtg_net = mtg_annotate_network(c)
