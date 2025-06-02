# %%
import os, datetime
import pandas as pd

# https://ndex2.readthedocs.io/
import ndex2.cx2

import yaml
import pickle
import pyshorteners

from microbetag.config import Config
from microbetag.networks import get_edgelist, read_cyjson

from microbetag.utils import (
    load_phenotypic_traits,
    extend_faprotax,
    extend_complements
)
from microbetag.seed_complementarity import load_seed_complement_files, build_url_with_seed_complements

# %% FUNCTIONS  -----------------------------------------------
def taxonomy_levels(node):
    # try:
    if len(node["v"]["microbetag::taxonomy"].split(";") ) == 7:
        (
            node["v"]["taxonomy::domain"],
            node["v"]["taxonomy::phylum"],
            node["v"]["taxonomy::class"],
            node["v"]["taxonomy::oder"],
            node["v"]["taxonomy::family"],
            node["v"]["taxonomy::genus"],
            node["v"]["taxonomy::species"]

        ) = node["v"]["microbetag::taxonomy"].split(";")


# %% Config  -----------------------------------------------

config_file = os.path.join(os.getcwd(),"test_data/test_build_cx2/config_v103.yml")
with open(config_file, 'r') as yaml_file:
    config = Config(yaml.safe_load(yaml_file), config_file)


# %% Co-occurrying taxa and their taxonomies  -----------------------------------------------

edgelist_df = get_edgelist(config)

# eg  {'nodeA': 'bin_185', 'nodeB': 'bin_192'}
taxa_pairs = edgelist_df[["nodeA", "nodeB"]].to_dict(orient="records")

# e.g.  'bin_31': 'd__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Reyranellales;f__Reyranellaceae;g__Reyranella;s__',
seq_id_to_taxonomy_dic = config.seq_to_taxon_df.set_index(
        config.seq_to_taxon_df.columns[0]
        )[ config.seq_to_taxon_df.columns[1] ].to_dict()



# %% Init attributes   -----------------------------------------------

nodes = []
node_names     = list(set(pd.concat([edgelist_df["nodeA"], edgelist_df["nodeB"]])))

for index, name in enumerate(node_names):
    node = {}
    node["id"] = index
    node["v"]  = {}
    node["v"]["name"]     = name
    node["v"]["microbetag::taxonomy"] = seq_id_to_taxonomy_dic[name]
    taxonomy_levels(node)
    nodes.append(node)

edges = []
for record in edgelist_df.iterrows():
    index, link = record
    node_a, node_b, weight = link
    edge = {}
    edge["id"] = index
    # NOTE (Haris Zafeiropoulos, 2025-03-26): It is the autoincreasing number to be used to map nodes to edges, not their names
    edge["s"]  = node_names.index(node_a)
    edge["t"]  = node_names.index(node_b)
    edge["v"]  = {}
    if weight > 0:
        edge["v"]["interaction type"] = "cooccurrence"
    else:
        edge["v"]["interaction type"] = "depletion"
    edge["v"]["microbetag::weight"]   = weight
    edges.append(edge)



# %% Phendb   -----------------------------------------------

# For each bin with a genome mapped, bin_phen_traits has a dictionary like:
# {'sulfate_reducer': {'presence': 'NO', 'confidence': 0.8646}, ..
# while, phentraits is just a set with all the available traits

def update_with_phen_traits(config, nodes):

    bin_phen_traits, _ = load_phenotypic_traits(config)

    for bin_name, phen_attributes in bin_phen_traits.items():
        # Get node based on the index of the node name
        node = nodes[node_names.index(bin_name)]

        if node["v"]["name"] == bin_name:
            for trait, values in phen_attributes.items():

                mtg_trait       = "::".join(["phendb", trait])
                mtg_trait_score = "::".join(["phendbScore", trait])

                node["v"][mtg_trait]       = values["presence"]
                node["v"][mtg_trait_score] = values["confidence"]
        else:
            print("This should not happening!! ")

update_with_phen_traits(config, nodes)

# %% FAPROTAX   -----------------------------------------------
def update_with_faprotax_traits(config, nodes):

    bin_faprotax_traits, _ = extend_faprotax(config)

    for bin_name, faprotax_attributes in bin_faprotax_traits.items():
        # Get node based on the index of the node name
        node = nodes[node_names.index(bin_name)]
        if node["v"]["name"] == bin_name:
            for trait in faprotax_attributes:
                mtg_trait            = "::".join(["faprotax", trait])
                node["v"][mtg_trait] = True
        else:
            print("This should not happening!! ")


update_with_faprotax_traits(config, nodes)

# %% MANTA   -----------------------------------------------

def update_with_manta(config, nodes):

    manta_net   = read_cyjson(config.manta_net)

    #  e.g. ('bin_179', 1.0),
    clusters = list(manta_net.nodes(data="cluster"))
    for cluster in clusters:
        bin_name, cluster = cluster
        node = nodes[node_names.index(bin_name)]
        node["manta::cluster"] = cluster

    # e.g.  ('bin_179', 'weak'),
    assignments = list(manta_net.nodes(data="assignment"))
    for assignment in assignments:
        bin_name, assignment = assignment
        node = nodes[node_names.index(bin_name)]
        node["manta::assignment"] = assignment

    # ('bin_184', {'x': -152.5582718971826, 'y': 237.3890227824764}),
    manta_layout = []
    positions    = list(manta_net.nodes(data="position"))
    for node_position in positions:
        bin_name, position = node_position
        manta_pos = {}
        manta_pos["node"]              = bin_name
        manta_pos["x"], manta_pos["y"] = position

    return manta_layout

# layout = update_with_manta(config, nodes)

# %%  Pathway Complementarities   -----------------------------------------------

def pathway_complement_edge(beneficiary, donor, complement):

    column     = ":".join(["compl:", beneficiary, donor])

    edge  = {}
    index = len(edges) + 1

    edge["id"] = index
    # NOTE (Haris Zafeiropoulos, 2025-03-26):
    # Conceptually, the donor is the source, since a compound would be secreted from it and drive to the beneficiary (target)
    edge["s"]  = node_names.index(donor)
    edge["t"]  = node_names.index(beneficiary)
    edge["v"]  = {}

    edge["v"]["interaction type"] = interaction_type
    edge["v"]["shared name"]      = " ".join([beneficiary, interacting, donor])
    edge["v"][column] = _hat_complement(complement)

    return edge


def _hat_complement(complements):

    hat_compl = []
    for compl in complements.values():
        hat_compl.append("^".join(compl))
    return hat_compl


interaction_type = "completes/competes with"
interacting      = "".join(["(", interaction_type, ")"])
complements_dict = extend_complements(
    complements_json=config.compl_file,
    descrps_path=config.module_descriptions,
    max_scratch_alt=config.max_scratch_alt,
    pathway_complements_dir=config.pathway_complements_dir,
    pathway_complement_percentage=config.pathway_complement_percentage,
)

nodes_in_compls_dict = set(complements_dict.keys())
for record in edgelist_df.iterrows():

    _, link           = record
    node_a, node_b, _ = link

    if node_a in nodes_in_compls_dict and node_b in nodes_in_compls_dict:

        # Node A the beneficiary - node B the donor
        complementAB = complements_dict[node_a][node_b]
        edges.append(pathway_complement_edge(node_a, node_b, complementAB))
        # Node B the beneficiary - node A the donor
        complementBA = complements_dict[node_b][node_a]
        edges.append(pathway_complement_edge(node_b, node_a, complementBA))


# %% Seed Complementarities  -----------------------------------------------


def verbose_seed_complement(complements, beneficiarys_nonseed, kmap, shortener):
    """
    Appends the seed complementarities between two taxa as attributes to their corresponding edge
    id_x:
    id_y:
    """

    maps_in         = list(kmap[kmap['modelseed'].isin(complements)]["map"].unique())
    complements_map = kmap[kmap['modelseed'].isin(complements)]
    beneficiarys_nonseeds_map = kmap[kmap['modelseed'].isin(beneficiarys_nonseed)]

    complements_verbose = []
    for kegg_map in maps_in:
        ksc = list(complements_map[complements_map["map"] == kegg_map]["kegg_compound"])
        msc = ";".join(set(complements_map[complements_map["map"] == kegg_map]["modelseed"]))
        ns = list(beneficiarys_nonseeds_map[beneficiarys_nonseeds_map["map"] == kegg_map]["kegg_compound"])
        surl = build_url_with_seed_complements(ksc, ns, kegg_map, shortener)
        des = kmap[kmap["map"] == kegg_map]["description"].unique().item()
        cat = kmap[kmap["map"] == kegg_map]["category"].unique().item()
        ksc = ";".join(set(ksc))
        complements_verbose.append([cat, des, msc, ksc, surl])

    merged_compl = ["^".join(gcompl) for gcompl in complements_verbose]

    return merged_compl





def seed_complement_edge(beneficiary, donor, complement, competition, cooperation, edges):

    # NOTE (Haris Zafeiropoulos, 2025-03-26):
    # Conceptually, the donor is the source, since a compound would be secreted from it and drive to the beneficiary (target)
    # This holds for the scores as well - for nodeA in scores, we consider its seeds. Thus, the A of the score should be the beneficiary, i.e. target

    source = node_names.index(donor)       # edge["s"]
    target = node_names.index(beneficiary) # edge["t"]
    column = f"seedCompl::{beneficiary}:{donor}"

    # Find existing edge
    edge = next((e for e in edges if e["s"] == source and e["t"] == target and e["v"]["interaction type"] == interaction_type), None)

    if edge is None:
        edge_index_to_replace = None
        index = len(edges) + 1
        edge = {
            "id": index,
            "s": source,
            "t": target,
            "v": {
                "interaction type": interaction_type,
                "shared name": f"{beneficiary} {interacting} {donor}"
            }
        }

    else:
        index = edge["id"]
        edge_index_to_replace = edges.index(edge)

    # Update edge attributes
    edge["v"].update({
        column: complement,
        "seed:competition": competition,
        "seed:cooperation": cooperation
    })

    return edge, edge_index_to_replace



shortener = pyshorteners.Shortener() if config.tinyurl else None

kmap = load_seed_complement_files(config.kegg_mappings)

seed_scores = pd.read_csv(config.phylomint_scores, sep="\t", header=None, skiprows=1)
seed_scores.columns = ["A", "B", "Competition", "Complementarity"]


with open(config.module_nonseeds, "rb") as f:
    non_seed_sets = pickle.load(f)

with open(config.seed_complements, "rb") as f:
    seed_complements = pickle.load(f)



seed_complements_dict = seed_complements.to_dict(orient="index")


nodes_in_seed_compls_dict = set(seed_complements_dict.keys())
for record in edgelist_df.iterrows():

    _, link           = record
    node_a, node_b, _ = link

    if node_a in nodes_in_seed_compls_dict and node_b in nodes_in_seed_compls_dict:

        # Node A the beneficiary - node B the donor
        scores = seed_scores.query('A == @node_a and B == @node_b')
        if not scores.empty:
            competAB, cooperAB = scores.iloc[0][["Competition", "Complementarity"]]
        else:
            competAB, cooperAB = None, None  # or suitable defaults

        seed_complementAB    = seed_complements_dict[node_a][node_b]
        beneficiarys_nonseed = non_seed_sets.loc[node_a].to_list()[0]

        seed_complementAB = verbose_seed_complement(seed_complementAB, beneficiarys_nonseed, kmap, shortener)

        e, update_edge_index =  seed_complement_edge(node_a, node_b, seed_complementAB, competAB, cooperAB, edges)
        if update_edge_index:
            edges[update_edge_index] = e
        else:
            edges.append(e)


        # Node B the beneficiary - node A the donor
        scores = seed_scores.query('A == @node_b and B == @node_a')
        if not scores.empty:
            competBA, cooperBA = scores.iloc[0][["Competition", "Complementarity"]]
        else:
            competBA, cooperBA = None, None  # or suitable defaults

        seed_complementBA  = seed_complements_dict[node_b][node_a]
        beneficiarys_nonseed = non_seed_sets.loc[node_b].to_list()[0]

        seed_complementBA = verbose_seed_complement(seed_complementBA, beneficiarys_nonseed, kmap, shortener)

        e, update_edge_index = seed_complement_edge(node_b, node_a, seed_complementBA, competBA, cooperBA, edges)
        if update_edge_index:
            edges[update_edge_index] = e
        else:
            edges.append(e)







# %% BUILD   -----------------------------------------------

# Create an empty net cx
net_cx = ndex2.cx2.CX2Network()

# Add nodes on the net_cx

for node in nodes:
    node_attributes = node["v"]
    # Add node
    net_cx.add_node(attributes=node_attributes)


# Add edges on the net_cx
for edge in edges:

    source = edge["s"]
    target = edge["t"]
    attributes = edge["v"].copy()

    if attributes["interaction type"] in ["depletion", "cooccurrence"]:
        attributes['microbetag::weight'] = float(edge["v"]['microbetag::weight'])

    filtered_attributes = {key: value for key, value in attributes.items() if not (isinstance(value, list) and len(value) == 0)}

    # create an edge connecting the nodes, id of edge is returned
    _ = net_cx.add_edge(source=source, target=target, attributes=filtered_attributes)



# if self.outfile is None:
    # Basename
timepoint = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M")
netfile = "_".join(["NEWmbtag_net", timepoint])
netfile = ".".join([netfile, "cx2"])
# if self.graphml_file is not None:
#     graphml_file = os.path.join(self.graphml_file, netfile)
# else:
graphml_file = netfile
# else:
#     graphml_file = self.outfile
net_cx.set_network_attributes({'name': 'microbetag annotated network'})
net_cx.write_as_raw_cx2(graphml_file)





# %%
# consider for layout..
# net_cx.set_visual_properties.__

"microbetag::ncbi-tax-level"

                {"po": edge_counter, "n": "shared name", "v": " ".join([id_a, "(cooccurss with)", id_b]), "d": "string"},
                {"po": edge_counter, "n": "interaction type", "v": "cooccurrence", "d": "string"}


                {"po": edge_counter, "n": "shared name", "v": " ".join([id_a, "(depletes)", id_b]), "d": "string"},
                {"po": edge_counter, "n": "interaction type", "v": "depletion", "d": "string"}


