"""
Aim:
    Annotate a microbial co-occurrence network with microbetag: phenotypic traits, pathway and seed
    complementarities, and network clusters.

Input (paths set in the configuration YAML):
    - The co-occurrence network as an edge list, and an abundance table or a sequence id to taxonomy map
    - The outputs of the external tools, as needed by the steps switched on: a KEGG ORTHOLOGY annotation
      file (pathway complementarity), GEMs (seed complementarity), phenotrex predictions, FAPROTAX
      sub-tables and a manta clustered network. The microbetag-nf pipeline produces all of them.

Output:
    - An annotated network in .cx2 format
"""

import os
import sys
import yaml
import argparse

from .utils import (
    mtg_logger,
    load_merged_ko_file,
)
from .config import load_config, Config
from .build_mtg_cx2 import mtg_annotate_network
from .seed_complementarity import run_seed_complementarity
from .pathway_complementarity import export_pathway_complementarities


logger = mtg_logger(__name__)


# This is the main function building a microbetag-annotated network
def run_microbetag(config: Config):
    """
    Main function for running the microbetag workflow.
    Based on the Config provided, it computes pathway and seed complementarities and annotates
    the network with them and with the outputs of the external tools.

    Returns:
        mtg_net: A microbetag-annotated network in CX2 format. CX2 is a JSON-based format, so it is
                easy to use for the response of the on-the-fly version to the query from MGG.
    Note:
        The mtg_net returned, is also saved as a .cx2 file in the output_directory using a timestamp
        on its filename, e.g. mtag_net_2025-05-08_17-47.cx2.
    """
    db = None
    if config.onthefly or config.api:
        from . import db
        from .helpers import otf_seqid_ncbi_gtdb_map

        # The onthefly config file should always include the credentials
        db.DB_CREDENTIALS = config.db_config

    # Phenotypic traits: from microbetagDB in the on-the-fly version; otherwise phenotrex predictions are an input
    _get_phen_annotations(config, db=db)

    # Maps required for otf in case of complementaritites
    if (config.path_compl or config.seed_compl) and config.onthefly:
        (
            config.pairs_of_interest,
            config.relative_genomes,
            config.mspecies_map_df

        ) = otf_seqid_ncbi_gtdb_map(config)

    # Pathway complementarity
    _get_path_compl(config, db=db)

    # Seed complementarity
    _get_seed_compl(config, db=db)

    # Annotate network in .cx format
    mtg_net = _annotate_network(config)

    return mtg_net


def _get_phen_annotations(config, db=None):
    """Phenotypic traits of the genomes from microbetagDB, in the on-the-fly version."""
    if config.phen_traits and config.onthefly:

        logger.info("⚙️ -- GETTING PHENOTYPIC TRAITS FROM MICROBETAGDB --")

        try:
            t = db.GetPhenotrexTraits(config)
            t.get_phen_traits()

        except Exception:
            error_msg = (
                "Phenotypic traits for the genomes under study failed to be exported from microbetagDB."
            )
            logger.error(error_msg)
            raise RuntimeError(error_msg)


def _get_path_compl(config, db=None):
    """Pathway complementarities, from the KEGG ORTHOLOGY annotation file or microbetagDB."""
    if not config.path_compl:
        return

    logger.info("⚙️ -- EXTRACTING PATHWAY COMPLEMENTARITIES --")

    if config.onthefly:
        db.get_path_compls_otf(config)

    elif not config.prev_path_compl:
        if not os.path.exists(config.alts_file) or not os.path.exists(config.pc_file):
            pivot_df = load_merged_ko_file(config.ko_merged)
            export_pathway_complementarities(config, pivot_df)


def _get_seed_compl(config, db=None):
    """Seed complementarities, from the GEMs or microbetagDB."""
    if not config.seed_compl:
        return

    logger.info("⚙️ -- EXTRACTING SEED COMPLEMENTARITIES --")

    if config.onthefly:

        logger.info("[INTERMEDIATE STEP] GET SEED COMPLS FROM MICROBETAG DB")

        config.get_scores      = True
        config.get_complements = True

        # Get dictionary with GTDB accession ids to their correspoding PATRIC
        gc_to_patric_ids        = db.patric_from_gc_list(config.repr_genomes_present)
        config.gc_to_patric_ids = db.update_for_patric(config.module_nonseeds, gc_to_patric_ids)

    run_seed_complementarity(config)


# Functions to support the main()
def _print_help():
    help_message = """
    Usage: microbetag --config <path_to_config_yml>

    Get the configuration template of this version with:
      microbetag --template config.yml

    Other options:
    -h        Display this help message.
    -v        Display version.
    """
    print(help_message)


def _print_version():

    from . import __version__

    print(f"microbetag version: {__version__}")


def _print_config_message():
    """Error message for failure during parsing the configuration YAML file."""

    conf_message = """
    The config file you provided cannot be parsed properly.
    Please make sure you follow the instructions on the documentation site:
    https://microbetag.readthedocs.io/en/latest/tutorials_local/config.html
    Start from the template of the microbetag version you are running: microbetag --template config.yml
    """
    logger.error(conf_message)


def _annotate_network(config):
    """Helper function to build annotated netowork based on microbetag findings."""
    mtg_net = None  # define upfront
    if config.precalc_only is False:
        logger.info("⚙️ -- ANNOTATE NETWORK ")
        mtg_net = mtg_annotate_network(config)

    # Keep arguments
    config.export_to_log(os.path.join(config.outdir, "parameters.log"))
    logger.info("A parameters.log file with the parameters used in this run was built.")

    logger.info("microbetag completed.")

    return mtg_net

# The main
def write_template(path: str) -> None:
    """Writes the configuration template of the installed microbetag version to `path`."""
    from importlib import resources

    if os.path.exists(path):
        sys.exit(f"{path} already exists; choose another file name.")

    template = resources.files("microbetag").joinpath("templates", "config.yml").read_text()
    with open(path, "w") as f:
        f.write(template)
    print(f"Configuration template written to {path}")


def main(argv=None):
    """
    Loads and parses a configuration YAML file 
    and invokes the main function for running the microbetag pipeline.
    """
    parser = argparse.ArgumentParser(description="Microbetag CLI")

    parser.add_argument("--config", "-c", help="Path to the configuration yaml file.")
    parser.add_argument(
        "--template", metavar="FILE",
        help="Write the configuration template of this microbetag version to FILE, then fill it in.",
    )
    parser.add_argument(
        "-v", "--version", action="store_true", help="Show Microbetag version"
    )

    args = parser.parse_args(argv)

    if args.version:
        _print_version()
        sys.exit()

    if args.template:
        write_template(args.template)
        sys.exit()

    elif args.config is None:
        _print_help()
        sys.exit(0)

    try:
        config = Config(load_config(args.config))

    except yaml.YAMLError:
        _print_config_message()
        sys.exit(1)

    # Run microbetag pipeline
    run_microbetag(config=config)


if __name__ == "__main__":

    main()
