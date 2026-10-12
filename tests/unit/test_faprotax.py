"""FAPROTAX on an abundance table: the merged functional table and one sub-table per trait."""

import os

import microbetag
from microbetag.helpers import Faprotax
from microbetag.tools import run_faprotax


class Config:
    def __init__(self, abundance_table, outdir):
        # Package directory, where mtg_maps_models/FAPROTAX_* lives
        self.cwd                  = os.path.dirname(microbetag.__file__)
        self.outdir               = str(outdir)
        self.abundance_table      = str(abundance_table)
        self.taxonomy_column_name = "taxonomy"
        self.delimiter            = "\t"
        self.__dict__.update(vars(Faprotax(config=self)))


def test_run_faprotax(tmp_path, data_dir):
    config = Config(data_dir / "faprotax" / "input_files" / "thirty_Samples.tsv", tmp_path)

    run_faprotax(config=config)

    assert os.path.getsize(config.faprotax_funct_table) > 0
    assert any(f.endswith(".txt") for f in os.listdir(config.faprotax_sub_tables))
