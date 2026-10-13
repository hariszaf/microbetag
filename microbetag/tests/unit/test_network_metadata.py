"""Nodes of metadata variables: the metadata_file of the configuration, and the metavar mark of their nodes."""


import pandas as pd
import yaml

from microbetag import Config
from microbetag.build_mtg_cx2 import init_nodes_and_edges

METAVAR = ["metavar"]
LEVEL = "microbetag::ncbi-tax-level"


def write_metadata(path):
    # One variable per row, its name first; ENV4 is categorical (FlashWeave makes a node per value: ENV4_A)
    path.write_text("ENV1\t0\t1\t0\nENV4\tA\tB\tA\n")
    return path


def test_metavar_nodes(tmp_path):
    metadata = write_metadata(tmp_path / "metadata.tsv")
    edgelist = pd.DataFrame(
        [["bin_101", "bin_151", 0.6], ["bin_101", "ENV1", 0.3], ["bin_151", "ENV4_A", -0.2]],
        columns=["node_A", "node_B", "weight"],
    )
    taxonomy = {"bin_101": "Bacteria;Firmicutes", "bin_151": "Bacteria;Firmicutes"}

    nodes, names, edges = init_nodes_and_edges(edgelist, taxonomy, False, str(metadata))

    attrs = {n["v"]["name"]: n["v"] for n in nodes}
    assert attrs["ENV1"].get(LEVEL) == METAVAR
    assert attrs["ENV4_A"].get(LEVEL) == METAVAR
    assert attrs["bin_101"].get(LEVEL) != METAVAR
    assert len(edges) == 3


def test_config_reads_metadata_file(data_copy):
    work = data_copy("microbetag")
    write_metadata(work / "inputs" / "metadata.tsv")
    with open(work / "config.yml") as f:
        conf = yaml.safe_load(f)
    conf["metadata_file"] = {"file_path": "inputs/metadata.tsv"}
    with open(work / "config.yml", "w") as f:
        yaml.safe_dump(conf, f)

    config = Config(conf, str(work / "config.yml"))
    assert config.metadata_file == str(work / "inputs" / "metadata.tsv")
