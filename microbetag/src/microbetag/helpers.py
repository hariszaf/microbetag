"""
Handlers classes allowing the different steps
"""

import os
import sys
import json
import pandas as pd
from typing import TYPE_CHECKING

from .utils import (
    mtg_logger,
    safe_literal_eval,
    resolve_file_path,
    convert_to_json_serializable,
)
from .networks import build_base_graph, get_edgelist

if TYPE_CHECKING:
    from .config import Config

class Emojis:
    def __init__(self) -> None:
        self.WARNING_EMOJI     = "\u2757"
        self.TADA_EMOJI        = '"\U0001f389"'
        self.RED_CROSS_EMOJI   = "\u274c"
        self.GREEN_CHECK_EMOJI = "\u2705"
        self.ANNOUNCEMET       = "\u1f4e3"

_logger_ = mtg_logger(__name__)
emojis   = Emojis()


def missing_input_error(step: str, key: str, kind: str, produced_by: str) -> ValueError:
    """Error for a step that is switched on in the YAML, without the input it needs."""
    path_key = "dir_path" if kind == "folder" else "file_path"
    return ValueError(
        f"{step} is switched on, but `{key}: {path_key}` is not set in the configuration file.\n"
        f"microbetag uses the {kind} produced by {produced_by}; run that step with the microbetag-nf "
        f"pipeline (https://github.com/msysbio/microbetag-nf), or point `{key}` to its output."
    )


class PathwayComplementarity:
    """
    Sets variables regarding pathway complementarity tasks based on user's config (.yml) file

    Args:
        config: An instance of the :class:`.config.Config` class.

    Note:
        From the :class:`.config.Config` class, the following attributes are being used:
            - `base_dir`
            - `ouput_dir`
            - `prev_path_compl`
            - `complement_max_length`
            - `ko_merged_file`
            - `pc_percentage`
    """

    def __init__(self, config: "Config"):
        self.conf       = config
        self.base_dir   = config.base_dir
        self.outdir = config.outdir

        # Init
        self.initialize(config)

    def initialize(self, conf):
        """Main function to initialize pathway complementarity settings."""

        # Get pathways
        self.output_dirs()

        # 3-column file with the KEGG ORTHOLOGY terms of each genome
        self.ko_merged = None
        ko_merged = (self.conf.yaml.get("ko_merged_file") or {}).get("file_path")

        if ko_merged:
            self.ko_merged = resolve_file_path(self.base_dir, ko_merged)

        elif not conf.onthefly and not self.prev_path_compl:
            raise missing_input_error(
                "Pathway complementarity", "ko_merged_file", "KEGG ORTHOLOGY annotation file",
                "annotating the genomes with kofam (hmmsearch)",
            )

    def output_dirs(self):
        """Paths to output folders and files"""

        pc_file      = (self.conf.yaml.get("pc_file") or {}).get("file_path")
        self.pc_dir     = os.path.join(self.outdir, "pathway_complementarity")
        self.pc_percent = (self.conf.yaml.get("pc_percent") or {}).get("value", 1)

        os.makedirs(self.pc_dir, exist_ok=True)

        if pc_file:

            self.pc_file      = resolve_file_path(self.base_dir, pc_file)
            self.prev_path_compl = True

            _logger_.info(f">>compl : {self.pc_file}")

        else:

            self.prev_path_compl = False
            self.pc_file      = os.path.join(self.pc_dir, "pathCompls.json")
            self.alts_file       = os.path.join(self.pc_dir, "alts.json")


class MappingPaths:
    """
    Sets paths to mapping files
    """

    def __init__(self):
        mtg                     = os.path.dirname(__file__)
        kegg_mappings           = os.path.join(mtg, "mtg_maps_models/kegg_mappings/")
        self.kegg_mappings      = kegg_mappings
        self.metanetx_compounds = os.path.join(mtg, "mtg_maps_models/MetaNetX/chem_xref.tar.gz")

        self.ref_ko_per_module = os.path.join(
            kegg_mappings, "kegg_terms_per_module.tsv"
        )
        self.modules_definitions_json_map = os.path.join(
            kegg_mappings, "module_definition_map.json"
        )
        self.kegg_modules_to_maps = os.path.join(kegg_mappings, "module_map_pairs.tsv")
        self.seed_ko_mo           = os.path.join(self.kegg_mappings, "seedId_keggId_module.tsv")
        self.module_descriptions  = os.path.join(kegg_mappings, "module_descriptions")


class Faprotax:
    """Path to the FAPROTAX sub-tables (one per trait) used to annotate the network nodes."""

    def __init__(self, config):
        sub_tables = (config.yaml.get("faprotax_sub_tables") or {}).get("dir_path")
        if sub_tables is None:
            raise missing_input_error(
                "FAPROTAX annotation", "faprotax_sub_tables", "folder of sub-tables",
                "FAPROTAX (collapse_table.py -s)",
            )
        self.faprotax_sub_tables = resolve_file_path(config.base_dir, sub_tables)


class NetworkHandler:
    """
    Set network-related configuration variables based on whether a network is already available or not.
    Check whether all the sequence identifiers present in the network as nodes, are also among those
    of the abundance table.
    """

    def __init__(self, config):

        if config.network:
            self.process_network(config)
        elif not config.precalc_only:
            raise ValueError(
                "Please provide the co-occurrence network to annotate (`edge_list: file_path`). "
                "To infer one from your abundance table with FlashWeave, use the microbetag-nf pipeline "
                "(https://github.com/msysbio/microbetag-nf)."
            )

    def process_network(self, config):
        """Process network edgelist and check bin consistency."""

        f = pd.read_csv(config.network, sep="\t")

        bins_in_net = set(f.iloc[:, 0]).union(f.iloc[:, 1])  # Get unique bin names

        if config.abundance_table is None:
            self.seq_ids = (
                bins_in_net  # Store sequence IDs if no abundance table is provided
            )


class AbdTableHandler:
    """
    Handles processing and validation of the abundance table related variables.

    Args:
        config: Instance of the :class:`Config` class
    """

    def __init__(self, config):

        if config.abundance_table is not None:
            self.load_abundance_table(config)

    def load_abundance_table(self, config):
        """Reads and validates the abundance table."""
        df = pd.read_csv(config.abundance_table, sep=config.delimiter)

        # Identify last column as taxonomy column
        self.taxonomy_column_name = df.columns[-1]
        non_numeric = (
            pd.to_numeric(df[self.taxonomy_column_name], errors="coerce").isna().any()
        )

        if not non_numeric:
            _logger_.error(
                "Taxonomy is not provided in the abundance table; "
                "at least not in the last column of the file as expected."
            )
            sys.exit(0)

        # First column is assumed to contain sequence IDs (bins)
        self.sequence_id_column_name = df.columns[0]
        self.seq_ids                 = df.iloc[:, 0].tolist()


class SeedComplementarityHandler:
    """
    Handles variables related to the seed complementarity module.
    The genome-scale metabolic models (GEMs) are an input: a folder with one SBML (.xml) file per genome.

    Args:
        config: Instance of the :class:`Config` class
    """

    def __init__(self, config):

        self.base_dir   = config.base_dir
        self.seed_compl = True
        self.user_models = True

        self.seeds_paths(config)

        if config.onthefly:
            return

        models = (config.yaml.get("genome_scale_models") or {}).get("dir_path")
        if models is None:
            if self.skip_sets:
                # Scores and complements only need the previously computed sets
                self.genres = None
                return
            raise missing_input_error(
                "Seed complementarity", "genome_scale_models", "folder of GEMs (.xml)",
                "a reconstruction tool (CarveMe, gapseq or ModelSEEDpy)",
            )

        self.genres = resolve_file_path(self.base_dir, models)
        self._detect_namespace()

    def _detect_namespace(self):
        """
        Detects whether the GEMs use the BiGG (CarveMe, gapseq) or the ModelSEED namespace,
        from the first species of the first model.
        """
        from libsbml import readSBML

        model_files = sorted(f for f in os.listdir(self.genres) if f.endswith(".xml"))
        if not model_files:
            raise ValueError(f"No GEMs (.xml) found in {self.genres}.")

        model = readSBML(os.path.join(self.genres, model_files[0])).getModel()
        if model is None or model.getNumSpecies() == 0:
            raise ValueError(f"Could not read the metabolites of {model_files[0]}.")

        first_species = model.getSpecies(0).getId()
        if first_species.startswith("M_"):
            first_species = first_species[2:]

        if first_species.startswith("cpd"):
            self.namespace                 = "modelseed"
            self.genre_reconstruction_with = "modelseedpy"
        else:
            # BiGG compounds are mapped to ModelSEED for the complements
            self.namespace                 = "BiGG"
            self.genre_reconstruction_with = "carveme"

    def seeds_paths(self, config):

        # Directory for seeds complementarity
        self.seeds_outdir  = os.path.join(config.outdir, "seeds_complementarity")
        os.makedirs(self.seeds_outdir, exist_ok=True)

        # NOTE (Haris Zafeiropoulos, 2025-04-29): These 2 are supposed to be the .json files.
        self.prev_conf     = resolve_file_path(self.base_dir, (config.yaml.get("prev_conf") or {}).get("file_path"))
        self.prev_nonseeds = resolve_file_path(self.base_dir, (config.yaml.get("prev_nonseeds") or {}).get("file_path"))

        self.skip_sets     = self.prev_conf is not None and self.prev_nonseeds is not None

        self.sc_pkl      = os.path.join(self.seeds_outdir, "seed_complements.pkl")
        self.seed_scores = os.path.join(self.seeds_outdir, "seed_scores.tsv")

        # NOTE (Haris Zafeiropoulos, 2025-05-07):
        # If prev_conf and prev_nonseeds are None, then module_seeds and module_nonseeds will be built during the run
        # Otherwise, these neeed to
        self.module_seeds    = os.path.join(self.seeds_outdir, "kegg_module_related_seeds.pkl")

        if config.onthefly:
            self.module_nonseeds = (config.yaml.get("prev_nonseeds_module") or {}).get("file_path")
        else:
            self.module_nonseeds = os.path.join(
                self.seeds_outdir, "kegg_module_related_nonseeds.pkl"
            )


def manta_input_net(config: "Config") -> bool:
    """Wrapper for building the intermediate network CYJSON file; input for manta"""

    manta_input        = build_base_graph(config)
    manta_input_serial = convert_to_json_serializable(manta_input)
    with open(config.base_network_file, "w") as f:
        json.dump(manta_input_serial, f, indent=4)
    return True


def otf_seqid_ncbi_gtdb_map(config: "Config") -> tuple[dict, dict, pd.DataFrame]:
    """
    Builds a dataframe with sequence ids of nodes A and B found associated in the
    co-occurrence network followed by their corresponding NCBI Taxonomy ids and the
    representative GTDB genomes.

    Returns:
        A tuple containing:
            - pairs_of_interest: {("",""). ("","")}
            - relative_genomes: {ncbi_id: [gc, gc, gc], ..}
            - mspecies_map_df: A data frame with the sequence ids of the abundance data and
                their mapped NCBI Taxonomy IDs and their corresponding GTDB representative genomes.
    Note:
        Strictly for the on-the-fly version
    """
    edgelist_df = get_edgelist(config.network)
    seq_map_df  = config.otf_seq_tax_df.dropna(subset=["species_ncbi_id"]).copy()

    merged = edgelist_df.merge(seq_map_df[['microbetag_id', 'species_ncbi_id', 'gtdb_gen_repr']],
                               left_on='node_A', right_on='microbetag_id', how='left') \
        .rename(columns={'species_ncbi_id': 'species_ncbi_id_A',
                         'gtdb_gen_repr': 'gtdb_gen_repr_A'}) \
        .drop(columns='microbetag_id')

    # Merge again to get info for nodeB
    merged = merged.merge(seq_map_df[['microbetag_id', 'species_ncbi_id', 'gtdb_gen_repr']],
                          left_on='node_B', right_on='microbetag_id', how='left') \
        .rename(columns={'species_ncbi_id': 'species_ncbi_id_B',
                         'gtdb_gen_repr': 'gtdb_gen_repr_B'}) \
        .drop(columns='microbetag_id')

    # Clean and apply transformations
    filtered = merged.dropna(subset=['gtdb_gen_repr_A', 'gtdb_gen_repr_B']).copy()

    filtered[['gtdb_gen_repr_A', 'gtdb_gen_repr_B']] = filtered[
        ['gtdb_gen_repr_A', 'gtdb_gen_repr_B']].applymap(safe_literal_eval)

    filtered[["species_ncbi_id_A", "species_ncbi_id_B"]] = filtered[
        ["species_ncbi_id_A", "species_ncbi_id_B"]].astype(int)

    # Explode and drop duplicates
    exploded        = filtered.explode('gtdb_gen_repr_A').explode('gtdb_gen_repr_B').reset_index(drop=True)
    mspecies_map_df = exploded.drop_duplicates(subset=['gtdb_gen_repr_A', 'gtdb_gen_repr_B'])

    # Build pairs of interest and relative genomes
    pairs_of_interest = {
        (str(row['species_ncbi_id_A']), str(row['species_ncbi_id_B']))
        for _, row in mspecies_map_df.iterrows()
    }
    pairs_of_interest.update({(b, a) for a, b in pairs_of_interest})
    relative_genomes = {
        str(row['species_ncbi_id_A']): {str(row['gtdb_gen_repr_A'])}
        for _, row in mspecies_map_df.iterrows()
    }
    relative_genomes.update(
        {str(row['species_ncbi_id_B']): {str(row['gtdb_gen_repr_B'])}
         for _, row in mspecies_map_df.iterrows()}
    )

    outfile = os.path.join(config.outdir, "edge_map.tsv")
    mspecies_map_df.to_csv(outfile)

    return pairs_of_interest, relative_genomes, mspecies_map_df
