"""
Aim:
    Gets an edgelist as input along with a df mentioning the NCBI Tax ids and the GTDB ids of the corresponding taxa
    and based on the annotation types asked, it build a cx-format graph that is going to be the final return object of microbetag

Author:
    Haris Zafeiropoulos

Based on:
    https://github.com/msysbio/microbetagApp/blob/develop/services/web/microbetag/scripts/buid_cx_annotated_graph.py
"""
import os
import copy
import pickle
import json
import pandas as pd
from utils import read_cyjson, load_seed_complement_files, build_url_with_seed_complements, convert_to_json_serializable

def build_cx_annotated_graph(conf):
    """
    If network provided by the user, needs to be in a .tsv format with 3 columns, for example:
    node_a    node_b    weight
    bin_23    bin_21    0.785
    """
    # Load edgelist
    if conf.flashweave:
        edgelist = pd.read_csv(conf.network, sep="\t", skiprows=2)
    else:
        edgelist = pd.read_csv(conf.network, sep="\t")
    set_of_nodes = set(edgelist.iloc[:, :2].values.ravel())
    unique_associated_pairs = {tuple(row) for _, row in edgelist.iloc[:, :2].iterrows()}

    # Load abundance table
    abundance_table_df = pd.read_csv(conf.abundance_table, sep="\t")
    seq_id_to_taxonomy = abundance_table_df[[conf.sequence_column_name, conf.taxonomy_column_name]]
    seq_id_to_taxonomy_dic = dict(zip(seq_id_to_taxonomy.iloc[:,0], seq_id_to_taxonomy.iloc[:,1]))

    annotated_cx = []
    init = {}
    init["numberVerification"] = [{"longNumber": 281474976710655}]
    annotated_cx.append(init)

    metadata = {}
    metadata["metaData"] = [
        {"name": "cyTableColumn", "version": "1.0"},
        {"name": "nodes", "version": "1.0"},
        {"name": "edges", "version": "1.0"},
        {"name": "nodeAttributes", "version": "1.0"},
        {"name": "edgeAttributes", "version": "1.0"},
        {"name": "networkAttributes", "version": "1.0"},
        {"name": "cartesianLayout", "version": "1.0"},
    ]
    annotated_cx.append(metadata)

    # ===========
    # GET COLUMN NAMES FOR ALL TABLES
    # ===========
    table_columns = {}
    table_columns["cyTableColumn"] = []
    cyTableColumns = [
        # for node table mandatory
        {"applies_to": "node_table", "n": "@id", "d": "string"},
        {"applies_to": "node_table", "n": "name", "d": "string"},
        {"applies_to": "node_table", "n": "shared name", "d": "string"},
        {"applies_to": "node_table", "n": "display name", "d": "string"},
        {"applies_to": "node_table", "n": "microbetag::taxon name", "d": "string"},
        {"applies_to": "node_table", "n": "microbetag::namespace"},
        {"applies_to": "node_table", "n": "microbetag::taxonomy", "d": "string"},
        {"applies_to": "node_table", "n": "microbetag::gtdb-genomes", "d": "list_of_string"},
        # for edge table mandatory
        {"applies_to": "edge_table", "n": "shared name"},
        {"applies_to": "edge_table", "n": "name"},
        {"applies_to": "edge_table", "n": "interaction type"},
        {"applies_to": "edge_table", "n": "microbetag::weight", "d": "double"}  # flashweave score
    ]

    """Phen traits"""
    bin_phen_traits = {}
    prediction_files = [os.path.join(conf.predictions_path, file) for file in os.listdir(conf.predictions_path)]
    phentraits = set()
    for file in prediction_files:
        if os.path.getsize(file) == 0:
            continue
        trait = pd.read_csv(file, sep="\t", skiprows=1)
        trait_name = os.path.basename(file).split(".prediction.tsv")[0]
        trait_filtered = trait[trait['Trait present'].notna()]
        trait_dict = trait_filtered.to_dict(orient="records")
        for case in trait_dict:
            bin_id, _ = os.path.splitext(case["Identifier"])
            if bin_id not in bin_phen_traits:
                bin_phen_traits[bin_id] = {}
            phentraits.add(trait_name)
            bin_phen_traits[bin_id][trait_name] = {}
            bin_phen_traits[bin_id][trait_name]["presence"] = case["Trait present"]
            bin_phen_traits[bin_id][trait_name]["confidence"] = case["Confidence"]

    for term in phentraits:
        cyTableColumns.extend([
            {"applies_to": "node_table", "n": "::".join(["phendb", term]), "d": "boolean"},
            {"applies_to": "node_table", "n": "::".join(["phendbScore", "".join([term, "Score"])]), "d": "double"}
        ])

    """FAPROTAX traits"""
    bin_faprotax_traits = {}
    fapro_sub_tables = [os.path.join(conf.faprotax_sub_tables, file) for file in os.listdir(conf.faprotax_sub_tables)]
    for file in fapro_sub_tables:
        trait_name, _ = os.path.splitext(os.path.basename(file))
        trait = pd.read_csv(file, sep="\t", skiprows=1)
        bins_with_trait = trait[conf.sequence_column_name].dropna()
        for bin_id in bins_with_trait:
            bin_faprotax_traits.setdefault(bin_id, []).append(trait_name)

    faprotax_traits = list(flatten_list(bin_faprotax_traits.values()))
    for term in faprotax_traits:
        column = {"applies_to": "node_table", "n": "::".join(["faprotax", term]), "d": "boolean"}
        cyTableColumns.append(column)

    """PATHWAY complements"""
    complements_dict = json.load(open(conf.compl_file))
    descrps = pd.read_csv(conf.module_descriptions, sep="\t", header=None)
    descrps.columns = ["category", "moduleId", "description"]
    column_order = ["moduleId", "description", "category"]
    descrps = descrps[column_order]

    complements_dict_ext = copy.deepcopy(complements_dict)

    for beneficiary_bin, potential_donors in complements_dict.items():
        for potential_donor, compls in potential_donors.items():
            if compls:
                complements_dict_ext[beneficiary_bin][potential_donor] = {}
                for compl in compls:
                    module_id = compl[0][3:]
                    kos_to_get = compl[1]
                    complet_alt = compl[2]

                    if len(complet_alt) == len(kos_to_get) > conf.max_scratch_alt:
                        print("Following complement was removed from the annotations as too hard to happen based on user's settings:")
                        print(compl)
                        print("----")
                        continue

                    perce = len(kos_to_get) / len(complet_alt)
                    if perce > conf.pathway_complement_percentage:
                        continue

                    compl_str = [x if isinstance(x, str) else ";".join(x) for x in compl[1:]]
                    triplet = descrps[descrps["moduleId"] == module_id].values.tolist()[0]
                    complements_dict_ext[beneficiary_bin][potential_donor][len(complements_dict_ext[beneficiary_bin][potential_donor])] = triplet + compl_str

    with open("pathway_complements_extended.json","w") as f:
        json.dump(complements_dict_ext, f)

    for edge in unique_associated_pairs:
        edge_col = {"applies_to": "edge_table", "n": "".join(["compl::", edge[0], ":", edge[1]]), "d": "list_of_string"}
        cyTableColumns.append(edge_col)

    """SEED complements"""

    kmap = load_seed_complement_files(conf.kegg_mappings)
    with open(conf.module_related_non_seeds, "rb") as f:
        non_seed_sets = pickle.load(f)

    with open(conf.seed_complements, "rb") as f:
        seed_complements = pickle.load(f)
    seed_complements_dict = seed_complements.to_dict(orient="index")

    for beneficiary, donors in seed_complements_dict.items():
        for donor, complements in donors.items():
            if len(complements) > 0:
                cyTableColumns.append(
                    {"applies_to": "edge_table",
                        "n": "".join(["seedCompl::", beneficiary, ":", donor]),
                        "d": "list_of_string"}
                )

    """SEED scores"""
    cyTableColumns.extend([
        {"applies_to": "edge_table", "n": "::".join(["seed", "competition"]), "d": "double"},
        {"applies_to": "edge_table", "n": "::".join(["seed", "cooperation"]), "d": "double"}
    ])
    seed_scores =pd.read_csv(conf.phylomint_scores, sep="\t", header=None, skiprows=1)
    seed_scores.columns = ["A", "B", "Competition", "Complementarity"]

    """MANTA CLUSTERS"""
    if conf.network_clustering:
        cartesianLayout = {}; cartesianLayout["cartesianLayout"] = []
        m1 = {"applies_to": "node_table", "n": "::".join(["manta", "cluster"]), "d": "double"}
        m2 = {"applies_to": "node_table", "n": "::".join(["manta", "assignment"]), "d": "string"}
        cyTableColumns.append(m1)
        cyTableColumns.append(m2)
        manta_output_file = "/".join([out_dir, 'manta_annotated.cyjs'])
        manta_net = read_cyjson(manta_output_file)
        clusters = list(manta_net.nodes(data="cluster"))
        assignments = list(manta_net.nodes(data="assignment"))
        positions = list(manta_net.nodes(data="position"))
        manta_annotations = {}
        for pair in clusters:
            manta_annotations[pair[0]] = {}
            manta_annotations[pair[0]]["cluster"] = pair[1]
        for pair in assignments:
            manta_annotations[pair[0]]["assignment"] = pair[1]
        for pair in positions:
            manta_annotations[pair[0]]["position"] = pair[1]

    table_columns["cyTableColumn"] = cyTableColumns
    annotated_cx.append(table_columns)

    # ===========
    # NETWORK TABLE
    # ===========
    networkAttributes = {}
    networkAttributes["networkAttributes"] = []
    networkAttributes["networkAttributes"].extend([
        {"n": "network type", "v": "microbetagAnnotated"},
        {"n": "name", "v": "microbetagNetwork"},
        {"n": "uri", "v": "https://hariszaf.github.io/microbetag/"},
        {"n": "version", "v": "1.0"}
    ])
    annotated_cx.append(networkAttributes)

    # NODES TABLE
    nodes = {}; nodes["nodes"] = []
    nodeAttributes = {}; nodeAttributes["nodeAttributes"] = []
    node_counter = 1000
    seq_to_nodeID = {}

    # ===========
    # NODES TABLE
    # ===========
    for bin_seq in set_of_nodes:

        node_counter += 1
        node = {"@id": node_counter, "n": bin_seq}
        seq_to_nodeID[bin_seq] = node_counter
        nodes["nodes"].append(node)

        taxonomy = seq_id_to_taxonomy_dic[bin_seq]

        # Node attributes
        nodeAttributes["nodeAttributes"].append({"po": node_counter, "n": "@id", "v": bin_seq, "d": "string"})
        nodeAttributes["nodeAttributes"].append({"po": node_counter, "n": "shared name", "v": bin_seq, "d": "string"})
        nodeAttributes["nodeAttributes"].append({"po": node_counter, "n": "display name", "v": bin_seq, "d": "string"})
        nodeAttributes["nodeAttributes"].append({"po": node_counter, "n": "microbetag::taxon name", "v": taxonomy.split(";")[-1], "d": "string"})
        nodeAttributes["nodeAttributes"].append({"po": node_counter, "n": "microbetag::taxonomy", "v": taxonomy, "d": "string"})

        # # Clusters
        # if conf.network_clustering:
        #     if bin_seq in manta_annotations.keys():
        #         nodeAttributes["nodeAttributes"].append({"po": node_counter, "n": "manta::cluster", "v": str(manta_annotations[seq]["cluster"]), "d": "double"})
        #         nodeAttributes["nodeAttributes"].append({"po": node_counter, "n": "manta::assignment", "v": manta_annotations[seq]["assignment"], "d": "string"})
        #         cartesianLayout["cartesianLayout"].append({"node": node_counter, "x": manta_annotations[seq]["position"]["x"], "y": manta_annotations[seq]["position"]["y"]})

        # phen traits
        if bin_seq in bin_phen_traits:
            for trait, presence in bin_phen_traits[bin_seq].items():
                trait_presence = True if presence["presence"] == "YES" else False
                nodeAttributes["nodeAttributes"].extend([
                    {"po": node_counter, "n": "::".join(["phendb", trait]), "v": trait_presence, "d": "boolean"},
                    {"po": node_counter, "n": "::".join(["phendbScore", "".join([trait, "Score"])]), "v": str(presence["confidence"]), "d": "double"}
                ])

        # FAPROTAX
        if bin_seq in bin_faprotax_traits and len(bin_faprotax_traits[bin_seq]) > 0:
            for term in bin_faprotax_traits[bin_seq]:
                nodeAttributes["nodeAttributes"].append({"po": node_counter, "n": "::".join(["faprotax", term]), "v": True, "d": "boolean"})

    annotated_cx.append(nodes)
    annotated_cx.append(nodeAttributes)

    print(len(seq_to_nodeID))
    print(len(set_of_nodes))

    # ===========
    # EDGES TABLE
    # ===========
    edges = {}; edges["edges"] = []
    edgeAttributes = {}; edgeAttributes["edgeAttributes"] = []
    edge_counter = node_counter + 1000

    for case in edgelist.iterrows():
        id_a = case[1][0]
        id_b = case[1][1]
        net_id_a = seq_to_nodeID[id_a]  # get the bin_id and then get its corresponding id on the net
        net_id_b = seq_to_nodeID[id_b]
        score = case[1][2]
        edge = {"@id": edge_counter, "s": net_id_a, "t": net_id_b, "i": "cooccurrence/depletion"}
        edges["edges"].append(edge)

        edgeAttributes["edgeAttributes"].append({"po": edge_counter, "n": "microbetag::weight", "v": str(score), "d": "double"})
        if score > 0:
            edgeAttributes["edgeAttributes"].extend([
                {"po": edge_counter, "n": "shared name", "v": " ".join([id_a, "(cooccurss with)", id_b]), "d": "string"},
                {"po": edge_counter, "n": "interaction type", "v": "cooccurrence", "d": "string"}
            ])
        else:
            edgeAttributes["edgeAttributes"].extend([
                {"po": edge_counter, "n": "shared name", "v": " ".join([id_a, "(depletes)", id_b]), "d": "string"},
                {"po": edge_counter, "n": "interaction type", "v": "depletion", "d": "string"}
            ])

        edge_counter += 1

        """
        Edge for A -> B
        """
        pot_edge = {"@id": (edge_counter), "s": net_id_a, "t": net_id_b, "i": "comp_coop"}

        # Path complements A -> B
        add_edge_pathway_complements(id_a, id_b, complements_dict_ext, edges, edgeAttributes, pot_edge, edge_counter)

        # Seed complements A -> B
        add_edge_seed_complements(id_a, id_b, seed_complements_dict, edges, edgeAttributes, pot_edge, edge_counter, kmap, non_seed_sets, seed_complements)

        # Seed scores A -> B
        add_seed_edge_attributes(id_a, id_b, seed_scores, edgeAttributes, edge_counter)

        """
        Edge for B -> A
        """
        pot_edge = {"@id": (edge_counter), "s": net_id_b, "t": net_id_a, "i": "comp_coop"}

        # Path complements B -> A
        add_edge_pathway_complements(id_b, id_a, complements_dict_ext, edges, edgeAttributes, pot_edge, edge_counter)

        # Seed complements B -> A
        add_edge_seed_complements(id_b, id_a, seed_complements_dict, edges, edgeAttributes, pot_edge, edge_counter, kmap, non_seed_sets, seed_complements)

        # Seed scores B -> A
        add_seed_edge_attributes(id_b, id_a, seed_scores, edgeAttributes, edge_counter)

    annotated_cx.append(edges)
    annotated_cx.append(edgeAttributes)

    # if cfg["manta"]:
    #     annotated_cx.append(cartesianLayout)

    # POST-metadata
    post_metadata = {}
    post_metadata["metaData"] = []
    post_metadata["metaData"].extend([
        {"name": "nodeAttributes", "elementCount": len(nodeAttributes["nodeAttributes"]), "version": 1.0},
        {"name": "edgeAttributes", "elementCount": len(edgeAttributes["edgeAttributes"]), "version": 1.0},
        {"name": "cyTableColumn", "elementCount": len(table_columns["cyTableColumn"]), "version": 1.0},
        {"name": "edges", "elementCount": len(edges["edges"]), "idCounter": node_counter + 1000, "version": 1.0},
        {"name": "nodes", "elementCount": len(nodes["nodes"]), "idCounter": 1001, "version": 1.0},
        {"name": "networkPropernetworkAttributesties", "elementCount": len(networkAttributes["networkAttributes"]), "version": 1.0},
    ])

    # if cfg["manta"]:
    #     post_metadata["metaData"].append({"name": "cartesianLayout", "elementCount": len(cartesianLayout["cartesianLayout"]), "version": 1.0})

    annotated_cx.append(post_metadata)

    # Status
    status = {}; status["status"] = []
    status["status"].append({"error": "", "success": True})
    annotated_cx.append(status)

    return annotated_cx

def seqId_faprotax_functions_assignment(path_to_subtables):
    """
    Parse the sub tables of the faprotax analysis
    to assign the biological processes related to each sequence id
    """
    seqId_faprotax_assignments = {}
    for subtable_file in os.listdir(path_to_subtables):
        f = os.path.join(path_to_subtables, subtable_file)
        process_name = subtable_file.split(".")[0].replace("_", " ")
        table_file = open(f, "r")
        table_file = table_file.readlines()
        for line in table_file[2:]:
            seqId = line.split("\t")[1]
            if seqId not in seqId_faprotax_assignments:
                seqId_faprotax_assignments[seqId] = [process_name]
            else:
                seqId_faprotax_assignments[seqId].append(process_name)
    return seqId_faprotax_assignments


def flatten_list(lista, flat_list=[]):
    """
    Recursive function taking as input a nested list and returning a flatten one.
    E.g. ['GCF_003252755.1', 'GCF_900638025.1', 'GCF_003252725.1', 'GCF_003253005.1', 'GCF_003252795.1', ['GCF_000210895.1'], ['GCF_000191405.1']]
    becomes ['GCF_003252755.1', 'GCF_900638025.1', 'GCF_003252725.1', 'GCF_003253005.1', 'GCF_003252795.1'].
    """
    for i in lista:
        if isinstance(i, list):
            flatten_list(i, flat_list)
        else:
            flat_list.append(i)
    return set(flat_list)


def add_edge_pathway_complements(id_x, id_y, complements_dict_ext, edges, edgeAttributes, pot_edge, edge_counter):
    if id_x in complements_dict_ext and id_y in complements_dict_ext[id_x]:
        pathway_complements = complements_dict_ext[id_x][id_y]
        edges["edges"].append(pot_edge)
        edgeAttributes["edgeAttributes"].extend([
            {"po": edge_counter, "n": "shared name", "v": f"{id_x} (completes/competes with) {id_y}", "d": "string"},
            {"po": edge_counter, "n": "interaction type", "v": "completes/competes with", "d": "string"}
        ])
        attr = f"compl::{id_x}:{id_y}"
        merged_compl = ["^".join(gcompl) for gcompl in pathway_complements.values()]
        edgeAttributes["edgeAttributes"].append({"po": edge_counter, "n": attr, "v": merged_compl, "d": "list_of_string"})


def add_edge_seed_complements(id_x, id_y, complements_dict, edges, edgeAttributes, pot_edge, edge_counter, kmap, non_seed_sets, seed_complements):
    if id_x in complements_dict and id_y in complements_dict[id_x]:
        if pot_edge not in edges["edges"]:
            edges["edges"].append(pot_edge)
            edgeAttributes["edgeAttributes"].extend([
                {"po": edge_counter, "n": "shared name", "v": f"{id_x} (completes/competes with) {id_y}", "d": "string"},
                {"po": edge_counter, "n": "interaction type", "v": "completes/competes with", "d": "string"}
            ])

        complements = complements_dict[id_x][id_y]
        complements_map = kmap[kmap['modelseed'].isin(complements)]
        maps_in = list(kmap[kmap['modelseed'].isin(complements)]["map"].unique())

        beneficiarys_nonseed = non_seed_sets.loc[id_x].to_list()[0]
        beneficiarys_nonseeds_map = kmap[kmap['modelseed'].isin(beneficiarys_nonseed)]

        complements_verbose = []
        for kegg_map in maps_in:
            ksc = list(complements_map[complements_map["map"] == kegg_map]["kegg_compound"])
            msc = ";".join(set(complements_map[complements_map["map"] == kegg_map]["modelseed"]))
            ns = list(beneficiarys_nonseeds_map[beneficiarys_nonseeds_map["map"] == kegg_map]["kegg_compound"])
            surl = build_url_with_seed_complements(ksc, ns, kegg_map)
            des = kmap[kmap["map"] == kegg_map]["description"].unique().item()
            cat = kmap[kmap["map"] == kegg_map]["category"].unique().item()
            ksc = ";".join(set(ksc))
            complements_verbose.append([cat, des, msc, ksc, surl])

        attr = f"seedCompl::{id_x}:{id_y}"
        merged_compl = ["^".join(gcompl) for gcompl in complements_verbose]
        edgeAttributes["edgeAttributes"].append({"po": edge_counter, "n": attr, "v": merged_compl, "d": "list_of_string"})


def add_seed_edge_attributes(id_x, id_y, seed_scores, edgeAttributes, edge_counter):
    matching_rows = seed_scores[(seed_scores['A'] == id_x) & (seed_scores['B'] == id_y)]
    if not matching_rows.empty:
        comp = matching_rows["Competition"].item()
        coop = matching_rows["Complementarity"].item()
        edgeAttributes["edgeAttributes"].extend([
            {"po": edge_counter, "n": "seed::competition", "v": str(comp), "d": "double"},
            {"po": edge_counter, "n": "seed::cooperation", "v": str(coop), "d": "double"}
        ])


if __name__ == "__main__":
    """
    Build an annotated .cx using already build objects
    """
    import sys
    import yaml
    from config import Config
    config_file = sys.argv[1]

    with open(config_file, 'r') as yaml_file:
        config = Config(yaml.safe_load(yaml_file), config_file)

    annotated_network = build_cx_annotated_graph(config)

    with open("test.cx", "w") as f:
            annotated_network2file = convert_to_json_serializable(annotated_network)
            json.dump(annotated_network2file, f)


