"""
Command-line entry points used by the microbetag-nf Nextflow pipeline.

Each module exposes a ``main(argv=None)`` function registered as a console script in setup.py.
Positional arguments mirror the former ``nf/bin/mtg/*.py`` scripts so that pipeline modules
only need to swap the command name.
"""
import os

import microbetag

_MAPS_MODELS       = os.path.join(os.path.dirname(microbetag.__file__), "mtg_maps_models")
KEGG_MAPPINGS      = microbetag._KEGG_MAPPINGS
SEED_KEGG_MODULE   = os.path.join(KEGG_MAPPINGS, "seedId_keggId_module.tsv")
METANETX_COMPOUNDS = os.path.join(_MAPS_MODELS, "MetaNetX", "chem_xref.tar.gz")
