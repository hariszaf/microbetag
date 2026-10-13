import os
import ast

import numpy as np
import pandas as pd

from fuzzywuzzy import fuzz

from typing import TYPE_CHECKING

from functools import partial
from concurrent.futures import ProcessPoolExecutor

from microbetag import utils as mtg_utils
from microbetag import mtg_logger, build_mtg_cx2


if TYPE_CHECKING:
    from .app_utils import Mappings


logger = mtg_logger(__file__)


def map_seq_to_ncbi_tax_level_and_id(
    abd_tab: pd.DataFrame, tax_col: str, seqId: str,
    tax_scheme: str, tax_delim:str, mappings: "Mappings",
    get_children=None
) -> tuple[pd.DataFrame, list]:
    """

    Returns 
        splt_tax: A pd.DataFrame with taxonomy levels, NCBI accession ids and
                  GTDB representative genomes (if available). A file callsed seq_map.tsv
                  with its content is saved in the main output folder.
    Example:
        splt_tax:
        >>> df.iloc[4,:]
        Unnamed: 0                              5
        domain                           Bacteria
        phylum                         Firmicutes
        class                          Clostridia
        order                     Oscillospirales
        family                   Oscillospiraceae
        genus                      Flavonifractor
        species            Flavonifractor plautii
        seqId                             ASV0012
        microbetag_id                     ASV0012
        extendedSpecies    Flavonifractor plautii
        ref_species_name         Flavonifractor plautii
        species_ncbi_id                  411475.0
        gtdb_gen_repr         ["GCA_000239295.1"]
        genus_ncbi_id                    946234.0
        family_ncbi_id                   216572.0
        ncbi_tax_id                      411475.0
        ncbi_tax_level                   mspecies
        Name: 4, dtype: object
    """
    # To-do [2025-04-23, Haris Zafeiropoulos]:
    # Check if there is an easy way to go on with get_children

    splt_tax = split_and_validate_taxonomy(abd_tab, seqId, tax_col, tax_delim)

    splt_tax = normalize_species_column(splt_tax)

    splt_tax, gtdb_on = add_species_ncbi_ids(splt_tax, tax_scheme, mappings)

    splt_tax = add_higher_level_ncbi_id(
        splt_tax, "genus", mappings.GENERA_NCBI_IDS, underscore=True
    )
    splt_tax = add_higher_level_ncbi_id(
        splt_tax, "family", mappings.FAMILIES_NCBI_IDS, underscore=True
    )

    get_ncbi_tax_level(splt_tax)

    splt_tax.loc[:, "gtdb_gen_repr"] = splt_tax["gtdb_gen_repr"].apply(
        lambda x: str(x) if isinstance(x, list) else x
    )
    splt_tax = splt_tax.drop_duplicates()

    # Parse "gtdb_gen_repr" values safely
    splt_tax["gtdb_gen_repr"] = splt_tax["gtdb_gen_repr"].apply(safe_literal_eval)

    # Flatten and clean genome list
    repr_genomes_present = mtg_utils.remove_nan_from_list(
        mtg_utils.flatten_list(splt_tax["gtdb_gen_repr"])
    )

    return splt_tax, repr_genomes_present  # , None


def split_and_validate_taxonomy(abd_tab: pd.DataFrame, seqId: str, tax_col: str, tax_delim: str) -> pd.DataFrame:
    """
    microbetag will be usins a 7-level taxonomy in all cases.
    In case a taxonomy level is not available, microbetag will have it as None.
    """
    taxonomies = abd_tab[[seqId, tax_col, "microbetag_id"]]
    splt_tax   = taxonomies[tax_col].str.split(tax_delim, expand=True)

    if splt_tax[0].str.contains("Root").all():
        splt_tax = splt_tax.drop(columns=[0])

    if splt_tax.shape[1] != 7:
        raise ValueError("The taxonomy scheme provided is not a 7-level one.")

    splt_tax.columns = build_mtg_cx2.__RANKS__

    splt_tax[seqId]           = taxonomies[seqId].values
    splt_tax["microbetag_id"] = taxonomies["microbetag_id"].values

    return splt_tax


def normalize_species_column(splt_tax: pd.DataFrame) -> pd.DataFrame:
    """ """
    species = splt_tax["species"].apply(process_underscore_taxonomy)
    genus   = splt_tax["genus"].apply(process_underscore_taxonomy)

    if splt_tax["species"].str.contains("s__").all():

        if splt_tax["species"].str.contains(r"__[\s_]").all():
            return splt_tax.assign(extendedSpecies=species)

        if any(g in s for g, s in zip(genus, species)):
            return splt_tax.assign(extendedSpecies=species)

        return splt_tax.assign(extendedSpecies=genus + " " + species)

    return splt_tax.assign(extendedSpecies=splt_tax["species"])


def add_species_ncbi_ids(splt_tax: pd.DataFrame, tax_scheme: str, mappings: "Mappings") -> tuple[pd.DataFrame, bool]:
    """
    Map taxonomy provided to microbetag taxonomy scheme to map NCBI Taxonomy id and level.

    """

    if tax_scheme     == "GTDB":
        taxon_map_file = "gtdbSpecies2ncbiId2accession.tsv"

    elif tax_scheme   == "Silva":
        taxon_map_file = "gtdb_silvaSpecies2ncbi2accession.tsv"

    elif tax_scheme   == "microbetag_prep":
        taxon_map_file = "gc_accession_16s_gtdb_ncbid.tsv"

    else:
        taxon_map_file = "species2ncbiId2accession.tsv"

    # Load mapping file for taxonomy - NCBI Tax ID - assembly
    taxon_map_path      = os.path.join(mappings.TAXONOMY, taxon_map_file)
    tax_ncbi_gc         = pd.read_csv(taxon_map_path, sep="\t")
    tax_ncbi_gc.columns = ["ref_species_name", "species_ncbi_id", "gtdb_gen_repr"]

    # If taxonomy scheme to be used for exact match
    if tax_scheme in {"GTDB", "Silva", "microbetag_prep"}:

        splt_tax = pd.merge(
            splt_tax,
            tax_ncbi_gc,
            left_on="extendedSpecies",
            right_on="ref_species_name",
            how="left",
        )
        splt_tax["gtdb_gen_repr"] = splt_tax["gtdb_gen_repr"].apply(
            lambda x: [x] if pd.notna(x) else np.nan
        )
        return splt_tax, True

    else:

        # Fallback logic for fuzzy mapping if needed
        splt_tax = fuzzy_similarity(splt_tax, tax_ncbi_gc)

        return splt_tax, False


def add_higher_level_ncbi_id(splt_tax, level, mapping_file, underscore=True):
    """
    Will add column with NCBI Taxonomy Id of a taxonomic level higher than species/strain.
    """
    df = pd.read_csv(mapping_file, sep="\t", names=[level, "ncbi_tax_id"])
    if underscore:
        splt_tax[f"extended{level}"] = splt_tax[level].apply(
            process_underscore_taxonomy
        )
    else:
        splt_tax[f"extended{level}"] = splt_tax[level]

    splt_tax = pd.merge(
        splt_tax,
        df,
        left_on=f"extended{level}",
        right_on=level,
        how="left",
        suffixes=("", "_y"),
    ).rename(columns={"ncbi_tax_id": f"{level.lower()}_ncbi_id"})

    return splt_tax.drop(columns=[f"extended{level}", f"{level}_y"], errors="ignore")


def process_underscore_taxonomy(entry):
    """
    Convert "s__" like taxa to ncbi like ones
    """

    if entry is None:
        return ""

    if "__" not in entry:
        return entry.replace("_", " ")  # NOTE (Haris Zafeiropoulos, 2025-05-11): mirobetag_prep case
    parts = entry.split("__")
    if len(parts) == 1:
        return entry
    if len(parts[1]) > 0:
        return parts[1].replace("_", " ")
    else:
        return np.nan


def get_ncbi_tax_level(splt_tax):
    """
    Assign the lowest available taxonomic level for each sequence in a
    splt_tax DataFrame. The function will first check if the sequence has
    a species level assignment, then a genus level assignment, then a family
    level assignment. If none of the above are available, the ncbi_tax_id
    column will be left as NaN.

    Parameters
    ----------
    splt_tax : pd.DataFrame

    Returns
    -------
    pd.DataFrame
        The same DataFrame with the additional column ncbi_tax_id
        containing the lowest available taxonomic level assignment for
        each sequence. If no assignment could be made, the value will be
        NaN.
    """

    # Apply the taxonomic assignment function
    # Assign taxonomic ID using the assign_tax_id_for_node_level function
    splt_tax["ncbi_tax_id"] = splt_tax.apply(assign_tax_id_for_node_level, axis=1)

    # Convert ncbi_tax_id to numeric, handling errors by coercing them to NaN
    splt_tax["ncbi_tax_id"] = (
        pd.to_numeric(splt_tax["ncbi_tax_id"], errors="coerce")
        .astype("Int64")
        .astype(str)
    )

    splt_tax["ncbi_tax_level"] = splt_tax.apply(assign_tax_level, axis=1)

    return splt_tax


def assign_tax_id_for_node_level(row):
    # A generator expression that returns the first available non-NaN value
    return next(
        (
            row[level]
            for level in ["species_ncbi_id", "genus_ncbi_id", "family_ncbi_id"]
            if not pd.isna(row[level])
        ),
        np.nan,
    )


def assign_tax_level(row):
    # Check if "gtdb_gen_repr" is not NaN or empty, return "mspecies" if valid
    if pd.notna(row["gtdb_gen_repr"]):
        return "mspecies"  # Or any other relevant level you wish to assign

    # Proceed to check for the lowest available taxonomic level
    taxonomic_levels = ["species", "genus", "family"]
    taxonomic_columns = ["species_ncbi_id", "genus_ncbi_id", "family_ncbi_id"]

    # Iterate over the taxonomic levels and columns, returning the first non-NaN taxonomic level
    for level, column in zip(taxonomic_levels, taxonomic_columns):
        if not pd.isna(row[column]):
            return level

    return "higher than family"


def fuzzy_similarity(splt_tax, tax_ncbi_gc):
    """
    Gets a list of taxa names and checks what is their best hit against a ref taxonomy.
    """
    logger.info(splt_tax.head())
    logger.info(tax_ncbi_gc.head())

    unique_species  = splt_tax["species"].dropna().unique()
    mapping_records = tax_ncbi_gc.to_dict("records")

    # Run in parallel
    # NOTE (Haris Zafeiropoulos, 2025-05-16):
    # The partial() function provides a s way to "pre-fill" some arguments of a function
    # — like setting default values — so you can call it later with fewer parameters.
    find_best_match_partial = partial(find_best_match, records=mapping_records)

    with ProcessPoolExecutor() as executor:
        results = list(executor.map(find_best_match_partial, unique_species, chunksize=10))

    # Convert to DataFrame
    match_results = pd.DataFrame(
        results,
        columns=["query_species", "best_match", "species_ncbi_id", "similarity", "gtdb_gen_repr"]
    )

    high_conf_matches = match_results[match_results["similarity"] > 90]

    # Step 2: Merge with split_tax on the Species column
    split_tax_enriched = splt_tax.merge(
        high_conf_matches[
            ["query_species", "best_match", "species_ncbi_id", "gtdb_gen_repr"]
        ],
        how      = "left",
        left_on  = "species",
        right_on = "query_species"
    )

    # Optional: Drop "query_species" if not needed
    split_tax_enriched = split_tax_enriched.drop(columns=["query_species"])

    return split_tax_enriched


# Updated matching function
def find_best_match(species: str, records: dict) -> tuple[str, str, str, float, str]:
    """
    Uses the `fuzzywuzzy` library to get the closest NCBI Taxonomy Id and their corresponding
    GTDB representative genomes.

    Note:
        A threshold of 90 on `fuzz.ratio` score is then applied in the :func:`fuzzy_similarity`
        function, to only keep rather probable matches.

    Note:
        `records` columns names are expected to be:
        ```
        tax_ncbi_gc.columns = ["ref_species_name", "species_ncbi_id", "gtdb_gen_repr"]
        ```

    Attention:
        The `fuzzywuzzy` is now abandoned and replaced by `theFuzz` (https://github.com/seatgeek/thefuzz/).
        In the near future `microbetag` will update using `theFuzz`; its interface is a bit different.

    """
    buzz_taxa = {
        "uncultured bacterium"            : 77133,
        "uncultured organism"             : 155900,
        "uncultured beta proteobacterium" : 86027,
        "uncultured soil bacterium"       : 164851,
        "uncultured rumen bacterium"      : 136703,
        "uncultured gamma proteobacterium": 86473,
        "D_8__uncultured rumen protozoa"  : None,
        "uncultured delta proteobacterium": 34034,
    }

    if species in buzz_taxa:
        return (
            species,
            species,
            str(buzz_taxa[species]),
            100,
            None
        )

    # First try for exact match (case-insensitive)
    for record in records:
        if species.lower() == record["ref_species_name"].lower():
            return (
                species,
                record["ref_species_name"],
                record["species_ncbi_id"],
                100,
                record["gtdb_gen_repr"]
            )

    # If no exact match, do fuzzy matching
    best_record, score = max(
        ((rec, fuzz.ratio(species, rec["ref_species_name"])) for rec in records),
        key=lambda x: x[1]
    )

    return (
        species,
        best_record["ref_species_name"],
        best_record["species_ncbi_id"],
        score,
        best_record["gtdb_gen_repr"]
    )


def safe_literal_eval(x):
    if isinstance(x, str):
        x = x.strip()
        if x.startswith(("[", "{", "(", "'")):
            try:
                return ast.literal_eval(x)
            except (ValueError, SyntaxError):
                pass  # Fall through
    return x
