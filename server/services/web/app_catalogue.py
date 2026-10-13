"""
Taxonomy mapping for a server with its own genome catalogue (CATALOGUE=custom).

A custom microbetagDB is built from the precalc workflow of microbetag-nf and microbetag-export-db; its genomes
have no NCBI Taxonomy ids. Each taxon of the abundance table is matched to the closest genomes of the catalogue
(microbetag.taxonomy), from the genomeTaxonomy table, instead of going through NCBI Taxonomy and GTDB.

The returned table has the columns of map_seq_to_ncbi_tax_level_and_id (app_taxonomy.py), so that the rest of
the on-the-fly run is the same:

    - species_ncbi_id: a key per distinct set of matched genomes (not an NCBI Taxonomy id)
    - ncbi_tax_id: empty
    - ncbi_tax_level: the rank shared with the matched genomes (species, genus, ...)
    - gtdb_gen_repr: the matched genomes of the catalogue
"""

import os

import pandas as pd

from microbetag import db, mtg_logger
from microbetag.taxonomy import RANKS, GenomeTaxonomyIndex

logger = mtg_logger(__file__)

_INDEX = None


def custom_catalogue() -> bool:
    return os.environ.get("CATALOGUE", "gtdb").lower() == "custom"


def catalogue_index() -> GenomeTaxonomyIndex:
    """The genomes of the database and their taxonomies; loaded once per worker."""
    global _INDEX
    if _INDEX is None:
        rows = db.execute("SELECT genomeId, taxonomy FROM genomeTaxonomy;")
        _INDEX = GenomeTaxonomyIndex(dict(rows), min_rank=os.environ.get("CATALOGUE_MIN_RANK", "genus"))
        logger.info(f"Genome catalogue: {len(_INDEX)} genomes")
    return _INDEX


def map_seq_to_catalogue(
    abd_tab: pd.DataFrame, tax_col: str, seqId: str, tax_delim: str, index: GenomeTaxonomyIndex = None
) -> tuple[pd.DataFrame, list]:
    """
    Returns:
        The taxonomy table of the sequences (see the module docstring) and the matched genomes.
    """
    index = index or catalogue_index()

    rows, keys = [], {}
    for seq, mtg_id, lineage in abd_tab[[seqId, "microbetag_id", tax_col]].itertuples(index=False):
        ranks = [r.strip() for r in str(lineage).split(tax_delim)] if pd.notna(lineage) else []
        if ranks and ranks[0].lower() == "root":
            ranks = ranks[1:]
        ranks = (ranks + [None] * len(RANKS))[: len(RANKS)]

        match   = index.match(seq, ";".join(r for r in ranks if r))
        genomes = match.genomes or None
        key     = keys.setdefault(tuple(genomes), len(keys) + 1) if genomes else None

        rows.append([*ranks, seq, mtg_id, key, None, match.rank, genomes])

    df = pd.DataFrame(
        rows,
        columns=[*RANKS, seqId, "microbetag_id", "species_ncbi_id", "ncbi_tax_id", "ncbi_tax_level", "gtdb_gen_repr"],
    )
    matched = df["gtdb_gen_repr"].dropna()
    logger.info(f"{len(matched)} of {len(df)} taxa matched to genomes of the catalogue")

    return df, sorted({g for genomes in matched for g in genomes})
