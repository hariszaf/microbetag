"""
How the server matches the taxa of an abundance table to the genomes of its database (TAXON_MATCHING):

    species (default)  a public server, whose users come with any taxonomy (GTDB, SILVA, microbetag's prep): the
                       species name of each taxon is matched, exactly or nearly, to a GTDB species, through NCBI
                       Taxonomy (app_taxonomy.py); taxa without a species get no genomes
    lineage            a server of a user's own genomes, whose abundance table names taxa as its catalogue does:
                       each taxon is matched to the closest genomes of the database, those sharing its deepest
                       rank down to LINEAGE_MIN_RANK (genus by default), from the genomeTaxonomy table (this module)

The former variable, CATALOGUE (gtdb or custom) and CATALOGUE_MIN_RANK, still work.

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


MATCHINGS = {"species", "lineage"}
FORMER    = {"gtdb": "species", "custom": "lineage"}


def taxon_matching() -> str:
    """species or lineage (TAXON_MATCHING, or the former CATALOGUE); an unknown value is an error."""
    value = os.environ.get("TAXON_MATCHING", "").strip().lower()
    if not value and os.environ.get("CATALOGUE", "").strip():
        former = os.environ["CATALOGUE"].strip().lower()
        if former not in FORMER:
            raise ValueError(f"CATALOGUE={former}: use TAXON_MATCHING=species or TAXON_MATCHING=lineage")
        logger.warning(f"CATALOGUE={former} is TAXON_MATCHING={FORMER[former]} now")
        return FORMER[former]
    value = value or "species"
    if value not in MATCHINGS:
        raise ValueError(f"TAXON_MATCHING={value}: it is species (the default) or lineage")
    return value


def lineage_matching() -> bool:
    return taxon_matching() == "lineage"


def catalogue_index() -> GenomeTaxonomyIndex:
    """The genomes of the database and their taxonomies; loaded once per worker."""
    global _INDEX
    if _INDEX is None:
        rows = db.execute("SELECT genomeId, taxonomy FROM genomeTaxonomy;")
        min_rank = os.environ.get("LINEAGE_MIN_RANK") or os.environ.get("CATALOGUE_MIN_RANK") or "genus"
        _INDEX = GenomeTaxonomyIndex(dict(rows), min_rank=min_rank)
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
