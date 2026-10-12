"""Annotate a co-occurrence network with microbetag (microbetag-annotate)."""
import argparse
import sys
import traceback
from pathlib import Path

import pandas as pd
import yaml

from microbetag.cli import KEGG_MAPPINGS


class Config:
    def __init__(self, yaml_path: str):
        with open(yaml_path, "r") as f:
            data = yaml.safe_load(f)

        if not isinstance(data, dict):
            raise ValueError(f"YAML file {yaml_path} is not valid.")

        for key, value in data.items():
            setattr(self, key, value)

        self.kegg_mappings       = Path(KEGG_MAPPINGS)
        self.module_descriptions = self.kegg_mappings / "module_descriptions"

    def __repr__(self):
        attrs = ", ".join(f"{k}={v!r}" for k, v in self.__dict__.items())
        return f"{self.__class__.__name__}({attrs})"


def load_data(config):
    """Set the sequence id to taxonomy map (seq_to_taxon_df, seq_ids) on the config."""
    from microbetag.config import load_abundance
    from microbetag.networks import get_edgelist
    from microbetag.utils import detect_separator

    if config.network is None and config.abundance_file:
        (
            config.seq_to_taxon_df,
            config.sequence_id_column_name,
            config.taxonomy_column_name,
            config.delimiter,

        ) = load_abundance(config.abundance_file)

        config.seq_ids = config.seq_to_taxon_df["sequence_id"].unique().tolist()

    elif config.abundance_file is None and config.network:

        config.delimiter = detect_separator(config.sequence_taxonomy_map)

        seq_to_taxon_df         = pd.read_csv(config.sequence_taxonomy_map, sep=config.delimiter)
        seq_to_taxon_df.columns = ["sequence_id", "taxonomy"]
        config.seq_to_taxon_df  = seq_to_taxon_df
        config.seq_ids          = config.seq_to_taxon_df["sequence_id"].unique().tolist()

    elif config.abundance_file and config.network:

        # NOTE: Not all sequence ids in the seq_ids need to have a taxonomy in this case --
        # only those coming from the abundance table.
        network_df  = get_edgelist(config.network)
        net_seq_ids = (
            pd.concat([network_df.iloc[:, 0], network_df.iloc[:, 1]])
            .unique()
            .tolist()
        )

        (
            config.seq_to_taxon_df,
            config.sequence_id_column_name,
            config.taxonomy_column_name,
            config.delimiter,

        ) = load_abundance(config.abundance_file)

        abd_seq_ids    = config.seq_to_taxon_df["sequence_id"].unique().tolist()
        config.seq_ids = net_seq_ids + abd_seq_ids

    return config


def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(prog="microbetag-annotate", description=__doc__)
    parser.add_argument("--config_file", "-c", required=True, help="Path to YAML configuration file")
    parser.add_argument("--network", "-n", required=True, help="Path to network file")
    parser.add_argument("--faprotax", "-f", help="Path to FAPROTAX subtables")
    parser.add_argument("--clustered", "-m", help="Path to manta clustered network")
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    from microbetag.build_mtg_cx2 import mtg_annotate_network

    try:
        c = Config(args.config_file)
        c.network = args.network

        if args.faprotax:
            c.faprotax_sub_tables = args.faprotax

        if args.clustered:
            c.manta_net = args.clustered

        load_data(c)
        mtg_annotate_network(c)

    except Exception as e:
        print(f"❌ Unexpected error: {e}", file=sys.stderr)
        traceback.print_exc()
        # Non-zero exit so that Nextflow marks the task as failed
        sys.exit(1)


if __name__ == "__main__":
    main()
