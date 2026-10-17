"""
Aim:
    Establishing connection and query on microbetagDB.
    Useful for the API and the on-the-fly version of microbetag, more specifically for the 
    phenotrex-based phenotypic traits and the pathway complementarity steps. 

Notes:
    Needs a hidden file called .env_dev.json with the database credentials.
"""

import json
import os
import pickle
import random
from typing import TYPE_CHECKING, Dict, Set, Tuple

import mysql.connector
import pandas as pd
from mysql.connector import pooling

from .utils import convert_to_json_serializable, mtg_logger

_logger_ = mtg_logger(__file__)
_KEGG_MAPPINGS_ = os.path.join(
    os.path.dirname(__file__), "mtg_maps_models", "kegg_mappings"
)

DB_CREDENTIALS = {}


# -----------
# Database
# -----------


def execute(phrase: str) -> list[tuple]:
    """
    Establish a database connection and perform an action
    """
    cnx    = mysql.connector.connect(**DB_CREDENTIALS)
    cursor = cnx.cursor()

    cursor.execute(phrase)

    rows = cursor.fetchall()

    cursor.close()
    cnx.commit()
    cnx.close()

    return rows


def execute_in_a_pool(cursor: mysql.connector, query: str):
    """
    Executes a query using a connection from the connection pool.
    """
    try:
        cursor.execute(query)
        result = [row for row in cursor]
        return result
    except mysql.connector.Error as err:
        print(f"Something went wrong: {err}")
        print(query)


def init_connection_pool() -> mysql.connector.pooling.MySQLConnectionPool:
    """
    Initiates a connection pool.
    A pool opens a number of connections and handles thread safety when providing connections to requesters.
    For more see: https://dev.mysql.com/doc/connector-python/en/connector-python-connection-pooling.html
    """
    connection_pool = pooling.MySQLConnectionPool(
        pool_name = "microbetagDB_pool",
        pool_size = 5,
        user      = DB_CREDENTIALS["user"],
        password  = DB_CREDENTIALS["password"],
        host      = DB_CREDENTIALS["host"],
        port      = int(DB_CREDENTIALS.get("port", 3306)),
        database  = DB_CREDENTIALS["database"],
    )
    return connection_pool


# --------
# Mapping
# --------
def gc_unify(gc_list):
    """
    Removes duplicates of a genome that has entries both as GCA and GCF in the db.
    """
    return list({x.replace("GCA", "GCF") for x in gc_list})


def alt_genome_prefix(gc):
    """
    Switches GCA prefic of a genome accession id to GCF and vice-versa.
    """
    if gc.startswith("GCA_"):
        gc_alt = gc.replace("GCA_", "GCF_")
    elif gc.startswith("GCF_"):
        gc_alt = gc.replace("GCF_", "GCA_")
    else:
        return 0
    return gc_alt


def get_genomes_for_ncbi_tax_id(ncbi_tax_id: int):
    """
    Get the genome IDs corresponding to a given NCBI Taxonomy ID from the microbetagDB.

    Example: 1281578
    """

    query = f"""
        SELECT genomeId
        FROM genome2taxNcbiId
        WHERE ncbiTaxId = {ncbi_tax_id};
    """.strip()

    genome_ids = execute(query)

    return {ncbi_tax_id: gc_unify([x[0] for x in genome_ids])}


def get_ncbi_tax_id_for_genome(gc_id: str):
    # "GCA_018819265.1"

    query = f"""
        SELECT ncbiTaxId
        FROM genome2taxNcbiId
        WHERE genomeId = '{gc_id}';
    """.strip()

    ncbi_tax_ids = execute(query)

    return {gc_id: list(ncbi_tax_ids)}


def patric_from_gc_list(gc_accession_list: list):
    """
    Gets a list of GC accession ids and returns a dictionary where the GC ids are the keys
    and their corresponding PATRIC ids are the values.

    ["GCA_003184265.1"]
    """
    gc_to_patric_dict = {}
    for gc in gc_accession_list:

        # GCA_/GCF_ accessions are looked up with both prefixes; the ids of a custom catalogue as they are
        gc_alt = alt_genome_prefix(gc) or gc

        query = f"""
            SELECT patricId
            FROM patricId2genomeId
            WHERE gtdbGenomeAccession = '{gc}'
            OR gtdbGenomeAccession = '{gc_alt}';
        """.strip()

        patricId = execute(query)

        # NOTE (Haris Zafeiropoulos, 2025-05-08): Some PATRIC ids have a different suffix in the metadata file
        # than in their corresponding assembies used for the seed complementarity step.
        # Yet, keeping the part before the dot (.) - using int(patricId[0][0]) cannot be a solution since in
        # many cases we have several strains for the same species, meaning, several after the dot parts.
        gc_to_patric_dict[gc] = patricId[0][0] if patricId and patricId[0] else None

    return gc_to_patric_dict


class GetPhenotrexTraits:

    def __init__(self, config=None):

        if config is not None:
            self.config = config

    def get_phen_traits(self):
        """
        Returns predictions for a list of genomes

        repr_genomes_present, config.phen_traits_dir
        """
        # Get predictions for each trait, for all genomes matched to the taxonomies given from the user.
        traits_per_genome = {
            genome: self.get_phendb_traits(genome)
            for genome in self.config.repr_genomes_present
        }

        # Build a df from the predictions dictionary
        df = pd.DataFrame(traits_per_genome)

        # Drop 'gtdbId' only if it exists in the index to avoid errors
        df = df.drop("gtdbId", errors="ignore")  # drop to make it safe if "gtdbId" isn't present.

        # Export predictions per trait to a 3-col file.
        GetPhenotrexTraits.to_csv(df, self.config.phen_traits_dir)

    def get_phendb_traits(self, gtdb_genome_id: str):
        """
        Get phenotypical traits based on phenDB classes based on its GTDB representative genome
        "GCA_018819265.1"
        """
        gtdb_id = gtdb_genome_id.strip()

        rows = execute(GetPhenotrexTraits.phen_query(gtdb_id))

        if len(rows) == 0:
            alt_gtdb_id = alt_genome_prefix(gtdb_genome_id)
            rows        = execute(GetPhenotrexTraits.phen_query(alt_gtdb_id))

        if len(rows) == 0:
            # _logger_.info(f"Genome {gtdb_genome_id} is not available in microbetagDB.")
            return 0

        query          = "SHOW COLUMNS FROM phenDB;"
        colnames       = [next(iter(x)) for x in execute(query)]
        genomes_traits = {i: j for i, j in zip(colnames, rows[0])}

        return genomes_traits  # gtdb_genome_id

    @staticmethod
    def phen_query(gtdb_id):
        """ Builds query for microbetagDB """

        query = f"""
            SELECT * FROM phenDB
            WHERE SUBSTRING_INDEX(gtdbId, '.', 1) = SUBSTRING_INDEX('{gtdb_id}', '.', 1);
        """.strip()

        return query

    @staticmethod
    def to_csv(df, outdir=None):
        """
        Saves phenotypic trait predictions of a trait for a set of genomes in a file, in a 3-column format.
        Identifier, Trait present, Confidence
        The values the 'Trait present' column may get, is 'YES, 'NO', and 'N/A'.
        """
        if outdir is None:
            outdir = os.getcwd()

        # Loop through traits (every 2 rows)
        for i in range(0, df.shape[0], 2):

            trait = df.index[i]

            presence, score = df.iloc[i], df.iloc[i + 1]

            output = pd.DataFrame({
                'Identifier'   : presence.index,
                'Trait present': presence.values,
                'Confidence'   : score.values
            })

            # Write to .tsv
            trait    = f"{trait}.prediction.tsv"

            with open(os.path.join(outdir, trait), "w") as f:

                f.write(f"# Trait: {trait}\n")

                f.write("Identifier\tTrait present\tConfidence\n")

                output.to_csv(f, sep="\t", index=False, header=False)


# --------
# Pathway complementarity
# --------

def seed_sets_of(gems) -> Tuple[dict, dict]:
    """
    The seed sets of GEMs from seedSets and nonSeedSets: {gem: {seed: confidence}} and {gem: [non-seeds]}.
    """
    gems = sorted({str(g) for g in gems if g})
    if not gems:
        return {}, {}
    confidence = {
        gem: {s.rsplit(":", 1)[0]: float(s.rsplit(":", 1)[1]) for s in seeds.split(",") if s}
        for gem, seeds in _select("SELECT gem, seeds FROM seedSets WHERE gem IN (%s);", gems)
    }
    nonseeds = {
        gem: [c for c in cpds.split(",") if c]
        for gem, cpds in _select("SELECT gem, nonSeeds FROM nonSeedSets WHERE gem IN (%s);", gems)
    }
    return confidence, nonseeds


def write_seed_files(gems, outdir: str) -> Tuple[str, str, str, Set[str]]:
    """
    Writes the seed sets of the GEMs of a network (seed_sets_of) as the files the seed complementarity step reads:
    confidences and non-seeds (.json.gz), and the KEGG module related non-seeds (.pckl, a pd.DataFrame with a
    list per GEM). Only these GEMs are loaded, whatever the size of the database.

    Returns:
        the three paths, and the GEMs found
    """
    import gzip
    from .helpers import MappingPaths
    from .seed_complementarity import get_kegg_module_related

    confidence, nonseeds = seed_sets_of(gems)
    found = set(confidence) & set(nonseeds)
    os.makedirs(outdir, exist_ok=True)
    conf_path     = os.path.join(outdir, "network_conf.json.gz")
    nonseeds_path = os.path.join(outdir, "network_nonseeds.json.gz")
    module_path   = os.path.join(outdir, "network_nonseeds_module_related.pckl")
    with gzip.open(conf_path, "wt") as f:
        json.dump({g: confidence[g] for g in found}, f)
    with gzip.open(nonseeds_path, "wt") as f:
        json.dump({g: nonseeds[g] for g in found}, f)
    module_cpds = get_kegg_module_related(MappingPaths().seed_ko_mo)
    module_nonseeds = pd.DataFrame.from_dict({g: [[c for c in nonseeds[g] if c in module_cpds]] for g in found})
    with open(module_path, "wb") as f:
        pickle.dump(module_nonseeds.T, f)
    _logger_.info(f"Seed sets of {len(found)} of {len(set(gems))} GEMs from the database")
    return conf_path, nonseeds_path, module_path, found


def get_path_compls_otf(config):
    """
    On the fly way to get pathway complementarities. Builds the pathCompls.json and
    the pathway_complements_extended.json files.

    The first file, in the stand-alone version lools like this:
    {"bin_101": {
        "bin_101": [],
        "bin_151": [["md:M00019", ["K00826"], ["K01652", "K01653", "K00053", "K01687", "K00826"],
                    "https://www.kegg.jp/kegg-bin/show_pathway?map00290/K00826%09%23EAD1DC/K01687%09%2300A898/"], ..]
        },..
    }
    where the elements regarding the complement and the alternative are actually lists, and not strings, while the latter:
    {"bin_101": {
        "bin_101": [],
        "bin_151": {"0": ["M00019", "Valine/isoleucine biosynthesis", "Branched-chain amino acid metabolism",
                            "K00826", "K01652;K01653;K00053;K01687;K00826", "https://...]
        }
    }

    To build those files for the on-the-fly version (otf), this function makes use of the config.otf_seq_tax_df
    attribute of the config, as returned by the taxonomy.py script.

    It creates a list of tuples with the NCBI Taxonomy Ids of found related taxa (in the edgelist)
    and using their corresponding GTDB representative genomes gets their complements.

    Arguments:
        config: Insance of the microbetat'g config class but as edited in the app.py script where the otf_seq_tax_df
                attribute is added

    Returns:
        mspecies_map_df: A pd.DataFrame with nodes names, NCBI Taxonomy ids and GTDB representative genomes for those
                        edges of the network that were both mapped to at least a GTDB genome.
    """

    # Ge the complements!
    compls = get_path_compls_for_ncbi_ids(config.relative_genomes, config.pairs_of_interest)

    # NOTE (Haris Zafeiropoulos, 2025-05-05):
    # Fix complements in the format set for the MGG
    # We will use the genome id here, even if the same file in the stand-alone has the sequence ids.
    compls_format = {}
    for j in compls.values():
        for k, v in j.items():
            ben_genome, don_genome = k
            if ben_genome not in compls_format:
                compls_format[ben_genome] = {}
            new_cases = []
            for case in v:                       # case is like [[...]]
                sub = case[0]
                formatted = [
                    sub[0],                      # module_id
                    sub[1].split(";"),           # essential_kos wrapped in list
                    sub[2].split(";"),           # all_kos wrapped in list
                    sub[3]                       # status
                ]
                new_cases.append(formatted)
            compls_format[ben_genome][don_genome] = new_cases

    # Serialize the formatted complements dictionary
    compls_serial = convert_to_json_serializable(compls_format)

    with open(config.pc_file, "w") as f:
        json.dump(compls_serial, f)


def get_path_compls_for_ncbi_ids(
    relative_genomes: dict[str, set[str]], pairs_of_interest: set[tuple[str, str]]
) -> dict:
    """
    The pathway complementarities of pairs of taxa, computed from the per-genome tables: for each pair of their
    genomes, the alternatives of the beneficiary's incomplete modules (genomeAlternatives, uniqueComplements) whose
    missing KOs the donor has (genomeModuleKos).

    Arguments:
        relative_genomes: A dictionary {"1260918": {"GCF_002102185.1"}, "1819566": {"GCF_009711525.1"}}
        pairs_of_interest: {("1260918", "1819566")}

    Returns:
        {(taxon_a, taxon_b): {(genome_a, genome_b): [[[module, complement KOs, KOs of the path, KEGG map url]], ..]}}
    """
    from .pathway_complementarity import pair_complements

    genomes = {g for a, b in pairs_of_interest for t in (a, b) for g in relative_genomes.get(t, [])}
    if not genomes:
        return {}

    _logger_.info(f"===> Pathway complementarities of {len(pairs_of_interest)} pairs of taxa, {len(genomes)} genomes")
    alternatives, module_kos = genome_pathway_rows(genomes)
    complements = unique_complements({i for ids in alternatives.values() for i in ids})

    result = {}
    for taxon_a, taxon_b in pairs_of_interest:
        for genome_a in relative_genomes.get(taxon_a, []):
            beneficiary = alternatives.get(genome_a)
            if not beneficiary:
                continue
            for genome_b in relative_genomes.get(taxon_b, []):
                donor = module_kos.get(genome_b)
                if not donor or genome_a == genome_b:
                    continue
                ids = pair_complements(beneficiary, donor, complements)
                if ids:
                    result.setdefault((taxon_a, taxon_b), {})[(genome_a, genome_b)] = [
                        [complements[i]["row"]] for i in ids
                    ]
    return result


def _genome_forms(genome: str) -> list:
    """A genome accession and its other prefix (GCA_ / GCF_), as either can be stored."""
    alt = alt_genome_prefix(genome)
    return [genome, alt] if alt else [genome]


def _select(query: str, values: list, chunk: int = 2000) -> list:
    """Run a SELECT ... IN (%s) query on chunks of values."""
    cnx_pool = init_connection_pool()
    connection = cnx_pool.get_connection()
    cursor = connection.cursor()
    rows = []
    try:
        for i in range(0, len(values), chunk):
            part = values[i: i + chunk]
            cursor.execute(query % ",".join(["%s"] * len(part)), part)
            rows.extend(cursor.fetchall())
    finally:
        cursor.close()
        connection.close()
    return rows


def genome_pathway_rows(genomes: set[str]) -> tuple[dict, dict]:
    """
    The rows of the genomes in genomeAlternatives (ids of the alternatives of their incomplete modules) and
    genomeModuleKos (KOs per module), keyed by the genome ids asked for.
    """
    stored = {form: genome for genome in genomes for form in _genome_forms(genome)}
    alternatives = {}
    for genome, ids in _select("SELECT genome, complementIds FROM genomeAlternatives WHERE genome IN (%s);",
                               list(stored)):
        alternatives[stored[genome]] = [int(i) for i in ids.split(",") if i]
    module_kos = {}
    for genome, module, kos in _select("SELECT genome, KoModuleId, kos FROM genomeModuleKos WHERE genome IN (%s);",
                                       list(stored)):
        module_kos.setdefault(stored[genome], {})[module] = set(kos.split(";"))
    return alternatives, module_kos


def unique_complements(ids: set[int]) -> dict:
    """
    The alternatives of uniqueComplements: {id: {"module", "missing" (set of KOs), "row": [module, missing KOs,
    KOs of the path, coloured KEGG map url]}}.
    """
    if not ids:
        return {}
    rows = _select("SELECT complementId, KoModuleId, complement, pathway FROM uniqueComplements "
                   "WHERE complementId IN (%s);", sorted(ids))
    coloured = build_kegg_urls([[[module, missing, path]] for _, module, missing, path in rows])
    return {
        int(cid): {"module": module, "missing": set(missing.split(";")), "row": c[0]}
        for (cid, module, missing, _), c in zip(rows, coloured)
    }


def build_kegg_urls(genome_pair_compls):
    """
    Takes as input the complements list between two genomes and
    build urls to colorify the related to the module kegg map based on the KO terms of the beneficiary (pink)
    and those it gets from the donor (green).

    Notes:
        Some modules do not belong to any map, e.g. https://www.kegg.jp/module/M00705. In these cases, we will have a N/A value in the url.
    """

    # NOTE (Haris Zafeiropoulos, 2025-05-05):
    # The stand-alone version has a similar function on the pathway_complementarity.py, consider making a consensus one.

    # Constants
    color_map_base_url   = "https://www.kegg.jp/kegg-bin/show_pathway?"
    present_kos_color    = "%09%23EAD1DC/"
    complement_kos_color = "%09%2300A898/"

    # Load module-to-map mappings
    module_map_file = os.path.join(_KEGG_MAPPINGS_, "module_map_pairs.tsv")
    df              = pd.read_csv(module_map_file, sep="\t", header=None, names=["module", "map"])
    df["module"]    = df["module"].str.replace("md:", "", regex=False).str.strip()
    df["map"]       = df["map"].str.strip()

    module_map = dict(zip(df["module"], df["map"]))
    updated_complements = []
    for [module_id, complement_str, ko_terms_str] in [list(c[0]) for c in genome_pair_compls]:
        beneficiary_kos    = []
        complement_kos     = []
        complement_kos_set = set(complement_str.split(";"))
        ko_terms           = ko_terms_str.split(";")

        for ko in ko_terms:
            if ko in complement_kos_set:
                complement_kos.append(ko + complement_kos_color)
            else:
                beneficiary_kos.append(ko + present_kos_color)

        # Construct URL or fallback
        url    = "N/A"
        map_id = module_map.get(module_id)
        if map_id:
            url = f"{color_map_base_url}{map_id}/" + "".join(beneficiary_kos + complement_kos)

        # Append new colored URL to the record
        updated_complements.append([[module_id, complement_str, ko_terms_str, url]])

    return updated_complements
