"""
The input network of manta (.cyjs), with node taxonomies taken either from a sequence-to-taxonomy
file or from an abundance table. manta itself runs in the microbetag-nf pipeline.
"""

import os

import pandas as pd
import pytest

from microbetag.config import load_abundance
from microbetag.helpers import manta_input_net


class NetConfig:
    """Network plus a sequence-to-taxonomy file."""

    def __init__(self, input_dir, outdir):
        seq_tax_map = pd.read_csv(input_dir / "seq2taxonomy.tsv", sep="\t", names=["sequence_id", "taxonomy"])

        self.outdir            = str(outdir)
        self.network           = str(input_dir / "edgelist.csv")
        self.seq_to_taxon_df   = seq_tax_map
        self.seq_ids           = seq_tax_map["sequence_id"].unique().tolist()
        # manta needs the .cyjs suffix
        self.base_network_file = str(outdir / "basenet.cyjs")


class AbdTableConfig:
    """Network plus an abundance table."""

    def __init__(self, input_dir, outdir):
        self.outdir          = str(outdir)
        self.abundance_table = str(input_dir / "thirty_Samples.tsv")
        self.network         = str(input_dir / "edgelist.csv")
        (
            self.seq_to_taxon_df,
            self.sequence_id_column_name,
            self.taxonomy_column_name,
            _,
        ) = load_abundance(self.abundance_table)
        self.seq_ids           = self.seq_to_taxon_df["sequence_id"].unique().tolist()
        self.base_network_file = str(outdir / "basenet.cyjs")


@pytest.fixture(params=["based_on_net", "based_on_abd_table"])
def config(request, tmp_path, data_dir):
    input_dir = data_dir / "manta" / "input_files" / request.param
    cls = NetConfig if request.param == "based_on_net" else AbdTableConfig
    return cls(input_dir, tmp_path)


def test_manta_input_net(config):
    assert manta_input_net(config)
    assert os.path.getsize(config.base_network_file) > 0
