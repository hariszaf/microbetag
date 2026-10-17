import tempfile

from flask import jsonify

from app_utils import ArgumentMissing, NotFoundError, MicrobetagVars, Mappings


from microbetag import mtg_logger, db, ExportSeedComplementarities


_logger_ = mtg_logger(__file__)
mtg_vars = MicrobetagVars()
mappings = Mappings()


# -------
# Mappings
# -------
def _stringify_gc_patric_ids(genomes: dict | str) -> dict:
    """
    Given a dictionary of genomes, return a new dictionary with the same keys,
    but with the values converted to strings, and only including the ones that
    are not None.
    """
    if isinstance(genomes, dict):
        genomes = list(genomes.values())[0]
    elif isinstance(genomes, str):
        genomes = list(genomes)
    else:
        raise TypeError(f"The following value of the 'genomes' variable is neither dict nor str. {genomes}")
    return {
        k: str(v)
        for k, v in db.patric_from_gc_list(genomes).items()
        if v is not None
    }


def _get_patric_ids(beneficiaryId: str, donorId: str, typeCategory: str) -> tuple[dict, dict]:
    """
    Given a beneficiaryId, donorId, and typeCategory, return a tuple of 2 dictionaries.
    The keys of the dictionaries are the IDs of the genomes, and the values are the corresponding PATRIC IDs.
    The typeCategory parameter can be either "ncbiTaxonomyIds" or "ncbiGenomeIds".
    If it is "ncbiTaxonomyIds", the function will get the genome IDs based on the NCBI taxonomy IDs.
    If it is "ncbiGenomeIds", the function will get the genome IDs directly from the input parameters.
    """
    if typeCategory == "ncbiTaxonomyIds":
        ben_genomes   = db.get_genomes_for_ncbi_tax_id(beneficiaryId)
        donor_genomes = db.get_genomes_for_ncbi_tax_id(donorId)

    elif typeCategory == "ncbiGenomeIds":
        ben_genomes   = beneficiaryId
        donor_genomes = donorId

    # Get PATRIC IDs
    # {'2184738': ['GCF_011366425.1', 'GCF_011364525.1', 'GCF_011375725.1', 'GCF_011375535.1']}
    # I make values() as a list, and get [0] to actually get the whole list

    patric_ben_ids   = _stringify_gc_patric_ids(ben_genomes)
    patric_donor_ids = _stringify_gc_patric_ids(donor_genomes)

    return patric_ben_ids, patric_donor_ids


def _post_seed_mappings(
    typeCategory    : str,
    beneficiaryId   : str,
    donorId         : str,
    gc_patric_ben   : str = None,
    gc_patric_donor : str | None = None,
    patric_ben_ids  : list | None = None,
    patric_donor_ids: list | None = None,
) -> dict:
    """
    Form with accession ids for the seed complements and scores API.
    """
    if typeCategory == "ncbiTaxonomyIds":
        return {
            "beneficiary_maps": {
                "ncbi_tax_id"      : beneficiaryId,
                "gc_and_patric_ids": gc_patric_ben,
            },
            "donor_maps": {
                "ncbi_tax_id"      : donorId,
                "gc_and_patric_ids": gc_patric_donor,
            },
        }
    if typeCategory == "ncbiGenomeIds":
        return {
            "beneficiary_maps": {
                "genome_accession": beneficiaryId,
                "patric_id"       : patric_ben_ids,
            },
            "donor_maps": {
                "genome_accession": donorId,
                "patric_id"       : patric_donor_ids,
            },
        }
    return {}


def _resolve_seed_ids(beneficiaryId: str, donorId: str, typeCategory) -> tuple[list[str], list[str], dict]:
    """
    Resolve the IDs of the given beneficiary and donor to a list of their corresponding PATRIC IDs
    and a dictionary of their mappings.

    Parameters
    ----------
    beneficiaryId : str
        The ID of the beneficiary genome.
    donorId : str
        The ID of the donor genome.
    typeCategory : str
        The type of the IDs. Can be either "ncbiTaxonomyIds", "ncbiGenomeIds", or "patricGenomeIds".

    Returns
    -------
    list, list, dict
        A tuple of two lists and a dictionary. The first list contains the PATRIC IDs of the
        beneficiary genome, the second list contains the PATRIC IDs of the donor genome, and the
        dictionary contains the mappings of the IDs.
    """
    if typeCategory == "patricGenomeIds":
        return [beneficiaryId], [donorId], None

    gc_patric_beneficiary, gc_patric_donor = _get_patric_ids(
        beneficiaryId, donorId, typeCategory
    )

    maps = _post_seed_mappings(
        typeCategory,
        beneficiaryId,
        donorId,
        gc_patric_beneficiary,
        gc_patric_donor,
        list(gc_patric_beneficiary.values()),
        list(gc_patric_donor.values()),
    )

    return (
        list(gc_patric_beneficiary.values()),
        list(gc_patric_donor.values()),
        maps
    )


def _post_path_compls(res: dict) -> dict:
    """
    Formatter function on how to display pathway complementarities on the API case.

    Arguments:
        res: {}

    Returns:
        response (dict): {}

    """
    response = {}
    outer_idx = 0
    for (b_tax, d_tax), genome_data in res.items():
        item = {
            "beneficiary-ncbi-tax-id": b_tax,
            "donor-ncbi-tax-id": d_tax,
            "complements": []
        }
        for (b_genome, d_genome), compls in genome_data.items():
            item["complements"].append({
                "beneficiary-genome": b_genome,
                "donor-genome": d_genome,
                "complements": compls
            })
        response[outer_idx] = item
        outer_idx += 1
    return response


def _validate_parameters(*params):
    """
    Validate that required parameters are provided.
    """
    for param, name in params:
        if not param:
            raise ArgumentMissing(f"Please provide {name}.")
    return None


def val_compl_res(retrieved, beneficiaryId, donorId, typeCategory, datatype):

    if len(retrieved) > 0:
        return jsonify(retrieved)
    else:
        raise NotFoundError(
            f"No {datatype} were found for beneficiary id '{beneficiaryId}'"
            f"and donor id '{donorId}' of category type '{typeCategory}'"
        )

# ------
# Pathway compls
# ------

def get_pathway_compls(beneficiaryId: str, donorId: str, typeCategory: str):
    """
    Exports the pathway complementarities between a beneficiary and a donor taxon.
    """
    _validate_parameters(
        (beneficiaryId, "beneficiaryId"),
        (donorId, "donorId"),
        (typeCategory, "typeCategory"),
    )

    if typeCategory == "ncbiTaxonomyIds":
        ben_genomes = set(next(iter(db.get_genomes_for_ncbi_tax_id(beneficiaryId).values())))
        donor_genomes       = set(next(iter(db.get_genomes_for_ncbi_tax_id(donorId).values())))
        ben_tax_id, donor_tax_id = beneficiaryId, donorId

    elif typeCategory == "ncbiGenomeIds":

        res                = db.get_ncbi_tax_id_for_genome(beneficiaryId)

        ben_tax_id = res.get(beneficiaryId, [[None]])[0][0]

        res = db.get_ncbi_tax_id_for_genome(donorId)

        donor_tax_id = res.get(donorId, [[None]])[0][0]

        ben_genomes, donor_genomes = {beneficiaryId}, {donorId}

    pairs_of_interest = {(ben_tax_id, donor_tax_id)}
    relative_genomes = {
        ben_tax_id: ben_genomes,
        donor_tax_id: donor_genomes
    }

    path_compls = db.get_path_compls_for_ncbi_ids(relative_genomes, pairs_of_interest)

    if not path_compls:
        raise LookupError("No pathway complementarities found.")

    return _post_path_compls(path_compls)


# ------
# Seed compls
# ------

class seedConfig:
    """
    For the API case (for now; maybe to be used for the onthefly but hard to)
    """

    def __init__(self, patric_ids, func):
        self.api           = True
        self.onthefly      = False
        self.patric_ids    = patric_ids
        self.seeds         = mtg_vars.UPLOAD_FOLDER  # outdir
        self.threads       = 2
        self.skip_sets     = True
        # The seed sets of these GEMs only, from the database
        outdir = tempfile.mkdtemp(dir=mtg_vars.UPLOAD_FOLDER, prefix="seeds_")
        self.prev_conf, self.prev_nonseeds, _, _ = db.write_seed_files(patric_ids, outdir)
        self.seed_ko_mo    = mappings.SEED_KO_MO

        if func == "scores":
            self.get_scores      = True
            self.get_complements = False

        elif func == "compls":
            self.get_scores      = False
            self.get_complements = True

        else:
            self.get_scores      = True
            self.get_complements = True


def all_seed_compls(beneficiaries, donors) -> dict:
    """
    Compute all pairwise seed complementarity scores between the beneficiary and donor organisms.

    Parameters
    ----------
    beneficiaries : list
        List of PATRIC IDs of the beneficiary organisms
    donors : list
        List of PATRIC IDs of the donor organisms

    Returns
    -------
    dict
        Dictionary of dictionaries, where the outer keys are the beneficiary PATRIC IDs,
        and the inner keys are the donor PATRIC IDs. The values are the seed complementarity scores.
    """

    all_compls = {}
    for beneficiary in beneficiaries:

        all_compls[beneficiary] = {}

        for donor in donors:

            c          = seedConfig(patric_ids=[beneficiary, donor], func="compls")
            c.api      = True
            c.onthefly = False
            cls = ExportSeedComplementarities(c)

            _, compls = cls.get_scores_and_compls()

            try:
                all_compls[beneficiary][donor] = compls[beneficiary][donor]
            except Exception:
                _logger_.info("No complements for this PATRIC ids combination")
                pass

    if len(all_compls) == 0:
        raise LookupError("No seed complements found.")

    return all_compls


def all_seed_scores(beneficiaries, donors) -> dict:
    """
    Compute all pairwise competition and cooperation scores between the beneficiary and donor organisms.

    Parameters
    ----------
    beneficiaries : list
        List of PATRIC IDs of the beneficiary organisms
    donors : list
        List of PATRIC IDs of the donor organisms
    conf : dict
        Dictionary of seed sets and their corresponding confidence scores
    nonseeds : dict
        Dictionary of non-seed sets

    Returns
    -------
    dict
        Dictionary of all pairwise competition and cooperation scores between the beneficiary and donor organisms
    """
    beneficiaries_and_donors = beneficiaries + donors

    c = seedConfig(patric_ids=beneficiaries_and_donors, func="scores")

    cls = ExportSeedComplementarities(c)

    scores, _  = cls.get_scores_and_compls()

    if scores is None:
        raise LookupError("No seed scores found.")

    all_scores = {"seed-scores": scores.to_dict(orient="index")}

    return all_scores


def seed_base_response(beneficiaryId, donorId, typeCategory, compute_fn) -> dict:
    """
    Processes beneficiary and donor genome IDs to compute seed complementarity data.

    Validates input parameters and resolves PATRIC IDs for the given IDs and type category.
    Loads required non-seed and confidence data, and then invokes the provided computation
    function to get results. If applicable, updates the result with mapping data.

    Args:
        beneficiaryId (str): Genome ID (NCBI Genome/Taxonomy ID or PATRIC ID) for the beneficiary.
        donorId (str): Genome ID (NCBI Genome/Taxonomy ID or PATRIC ID) for the donor.
        typeCategory (str): Category of ID ("ncbiGenomeIds", "ncbiTaxonomyIds", "patricGenomeIds").
                compute_fn (function): Function to compute results using the resolved PATRIC IDs
                (all_seed_compls, all_seed_scores).

    Returns:
        JSON response with computed results and optional mappings, or an error if validation fails.
    """

    _validate_parameters(
        (beneficiaryId, "beneficiaryId"),
        (donorId, "donorId"),
        (typeCategory, "typeCategory"),
    )

    patric_beneficiaries, patric_donors, maps = _resolve_seed_ids(
        beneficiaryId, donorId, typeCategory
    )

    _logger_.info(patric_beneficiaries)
    _logger_.info("---")
    _logger_.info(patric_donors)
    _logger_.info("---")
    _logger_.info(maps)

    result = compute_fn(patric_beneficiaries, patric_donors)  # , conf, nonseeds

    if typeCategory != "patricGenomeIds":
        result.update(maps)

    return result
