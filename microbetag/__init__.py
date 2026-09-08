import os

from .config import Config
from .helpers import (
    AbdTableHandler,
    BinsHandler,
    Faprotax,
    MappingPaths,
    NetworkHandler,
    PathwayComplementarity,
    SeedComplementarityHandler,
    manta_input_net,
    otf_seqid_ncbi_gtdb_map,
)
from .pathway_complementarity import (
    all_alternatives,
    all_complements,
    build_kegg_url,
    export_pathway_complementarities,
)
from .seed_complementarity import (
    ExportSeedComplementarities,
    build_url_with_seed_complements,
    kegg_module_related_intersect,
    load_seed_complement_files,
)
from .utils import (
    extend_complements,
    extend_faprotax,
    ko_list_parser,
    load_phenotypic_traits,
    merge_ko,
    mtg_logger,
)

try:
    from .db import (
        GetPhenotrexTraits,
        get_genomes_for_ncbi_tax_id,
        get_ncbi_tax_id_for_genome,
        get_path_compls_for_ncbi_ids,
        patric_from_gc_list,
    )
except Exception:
    print(
        "mysql-connector-python is not installed in the running environment."
        "Dependency and microbetag features required only for on-the-fly version."
    )


from .build_mtg_cx2 import mtg_annotate_network
from .genres import GEMSReconstruction
from .microbetag import run_microbetag
from .networks import build_base_graph, get_edgelist, read_cyjson
from .tools import (
    hmmsearch,
    kegg_annotation,
    phenotrex_genotype,
    phenotrex_predict,
    run_faprotax,
    run_flashweave,
    run_manta,
    run_prodigal,
    run_seed_complementarity,
)

_KEGG_MAPPINGS = os.path.join(
    os.path.dirname(__file__), "mtg_maps_models", "kegg_mappings"
)
_KEGG_TERMS_PER_MODULE = os.path.join(_KEGG_MAPPINGS, "kegg_terms_per_module.tsv")
_MODULE_DEFINITION_MAP = os.path.join(_KEGG_MAPPINGS, "module_definition_map.json")
_KEGG_MODULES_TO_MAPS = os.path.join(_KEGG_MAPPINGS, "module_map_pairs.tsv")

_MTG_PHEN_ENV = "mtg-phenotrex"
_MTG_MODELSEED_ENV = "mtg-modelseed"

__version__ = "1.0.4"
__license__ = "GNU GPL3"
__authors__ = ["Haris Zafeiropoulos <haris.zafeiropoulos@kuleuven.be>"]
__cite__ = (
    "Zafeiropoulos H, Michail Delopoulos EI, Erega A, Schneider A, Geirnaert A, Morris J, Faust K."
    "microbetag: simplifying microbial network interpretation through annotation, enrichment tests and metabolic complementarity analysis."
    "Genome Biol 26, 292 (2025)"
    "DOI: https://doi.org/10.1186/s13059-025-03769-2"
)
