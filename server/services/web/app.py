import os
import yaml
import traceback

from flask_caching import Cache
from flask import jsonify, request

from app_utils import (
    ArgumentMissing, NotFoundError, FileFormatError,
    Mappings, UserArgsVars, MicrobetagVars,
    init_app, with_traceback_on_errors,
    build_data_files, sync_yaml, args_and_abd_data,
    check1000, prep_abundance_table
)

from app_taxonomy import map_seq_to_ncbi_tax_level_and_id
from app_catalogue import lineage_matching, map_seq_to_catalogue, taxon_matching
from app_pipeline import pipeline_available, run_network_workflow

from api import (
    get_pathway_compls, val_compl_res,
    all_seed_compls, all_seed_scores, seed_base_response
)

from microbetag import Config, mtg_logger, GetPhenotrexTraits, db, run_microbetag

__version__ = "1.0.2"

# How taxa are matched to genomes (TAXON_MATCHING): checked at start, so that a wrong value stops the server
TAXON_MATCHING = taxon_matching()

# --------
# Fire app
# --------
# app, login = init_app()
app, login, jl = init_app()

db.DB_CREDENTIALS = login


cache  = Cache(app, config={"CACHE_TYPE": "simple", "CACHE_DEFAULT_TIMEOUT": 20})
logger = mtg_logger(__file__)


# --------
# Error handlers
# --------
@app.errorhandler(NotFoundError)
def handle_not_found_error(e):
    return jsonify({"error": str(e)}), 404

@app.errorhandler(ArgumentMissing)
def handle_missing_argumen_error(e):
    return jsonify({"error": str(e)}), 400

@app.errorhandler(FileFormatError)
def handle_file_format_error(e):
    return jsonify({"status": "error", "message": str(e)}), 400

@app.errorhandler(Exception)
def handle_all_errors(e):
    tb = traceback.format_exc()
    return jsonify({"error": str(e), "traceback": tb}), 500


# -----------------------------------------------------------------------------
# ------------------------            API             -------------------------
# -----------------------------------------------------------------------------

@app.route("/", methods=["GET", "POST"])
def index():
    message = """
    Hello friend!

    This is a web-server of the Lab of Microbial Systems Biology, providing some of our resources and tools.
    Currently:

        - <a href="https://microbetag.readthedocs.io/en/v1.0.3/" target="_blank">microbetag</a> and 
        - <a href="https://bacterial-growth.readthedocs.io" target="_blank">μGrowthDB</a>.

    Visit <a href="http://www.msysbiology.com/#tools" target="_blank">the lab website</a> for
    an overview of our research and tools.
    """
    return message.replace("\n", "<br>")


@app.route("/ncbiTaxId-to-genomeId", methods=["GET"])
def get_related_genomes():
    """
    API for getting NCBI assembly accession ids of GTDB representative genomes available on microbetagDB.

    Arguments:
        ncbiTaxId

    Example:
        http://localhost:1337/ncbiTaxId-to-genomeId?ncbiTaxId=537007
    """

    ncbiTaxId = request.args.get("ncbiTaxId")
    if not ncbiTaxId:
        raise ArgumentMissing("Please provide a ncbiTaxId.")

    try:

        q = db.get_genomes_for_ncbi_tax_id(ncbiTaxId)
        if len(q[ncbiTaxId]) == 0:
            raise NotFoundError(f"NCBI Taxonomy Id {ncbiTaxId} was not found in microbetagDB.")

        return jsonify(q)

    except RuntimeError as e:

        return jsonify({"error": str(e)}), 500


@app.route("/genomeId-to-ncbiTaxId", methods=["GET"])
def get_ncbi_tax_id():
    """
    API for getting NCBI assembly accession ids of GTDB representative genomes available on microbetagDB.

    Arguments:
        ncbiTaxId

    Example:
        http://localhost:1337/genomeId-to-ncbiTaxId?genomeId=GCF_002222595.2
    """

    genomeId = request.args.get("genomeId")

    if not genomeId:
        raise ArgumentMissing("Please provide a genomeId, e.g. GCA_011050685")

    try:

        q = db.get_ncbi_tax_id_for_genome(genomeId)

        if len(q[genomeId]) == 0:

            raise NotFoundError(f"NCBI Taxonomy Id {genomeId} was not found in microbetagDB.")

        return jsonify(q)

    except RuntimeError as e:

        return jsonify({"error": str(e)}), 500


@app.route("/phen-traits", methods=["GET"])
def phen_traits():
    """
    API for getting phenotypic traits of a GTDB representative genome through its NCBI assembly accession id.

    Arguments:
        ncbiAssemblyId

    Example:
        http://localhost:1337/phen-traits?ncbiAssemblyId=GCF_002222595.2

    """
    gtdb_genome_id = request.args.get("ncbiAssemblyId")

    if not gtdb_genome_id:
        raise ArgumentMissing("Please provide a ncbiAssemblyId.")

    # Get traits from microbetagDB
    pt     = GetPhenotrexTraits()
    traits = pt.get_phendb_traits(gtdb_genome_id)

    if not traits:
        message = f"""
        No Phen traits for the genome id {gtdb_genome_id}. \
        Make sure you are asking for a GTDB v202 representative genome.
        """
        raise NotFoundError(message)

    return jsonify(traits)


# -------------------------
# Pathway complementarities
# -------------------------

@app.route("/pathway-complements", methods=["GET"])
def genomes_complements():
    """
    Gets pathway complements that a donor can provide to a baneficiary genome using either NCBI Taxonomy Ids, or
    NCBI assembly accession ids, or PATRIC ids of the corresponding GTDB representative genomes.

    Arguments:
        beneficiaryId (str): Genome ID (NCBI Genome/Taxonomy ID or PATRIC ID)
        donorId       (str): Genome ID (NCBI Genome/Taxonomy ID or PATRIC ID)
        typeCategory  (str): One of "ncbiGenomeIds", "ncbiTaxonomyIds"

    Example:
        http://localhost:1337/pathway-complements?beneficiaryId=GCA_011050685.1&donorId=GCA_002289575.1&typeCategory=ncbiGenomeIds

    """
    beneficiaryId = request.args.get("beneficiaryId")
    donorId       = request.args.get("donorId")
    typeCategory  = request.args.get("typeCategory")

    try:
        compls = get_pathway_compls(beneficiaryId, donorId, typeCategory)
    except Exception as e:
        return jsonify({"error": str(e)}), 500

    return val_compl_res(
        compls, beneficiaryId, donorId, typeCategory, "pathway complementarities"
    )


# -------------------------
# Seed complementarities and scores
# -------------------------

@app.route("/seed-complements", methods=["GET"])
def seed_complements():
    """
    For a pair of NCBI Taxonomy or GTDB representative genomes IDs or PATRIC IDs
    returns the seed complementarities for all combinations between
    the beneficiary and the donor PATRIC ids found on microbetagDB.

    Arguments:

        beneficiaryId: Genome ID (NCBI Genome/Taxonomy ID or PATRIC ID)
        donorId      : Genome ID (NCBI Genome/Taxonomy ID or PATRIC ID)
        typeCategory : One of "ncbiGenomeIds", "ncbiTaxonomyIds", "patricGenomeIds"

    Examples:

        http://localhost:1337/seed-complements?beneficiaryId=182270&donorId=161920&typeCategory=ncbiTaxonomyIds

        http://localhost:1337/seed-complements?beneficiaryId=GCF_002980625.1&donorId=GCF_011375535.1&typeCategory=ncbiGenomeIds

        http://localhost:1337/seed-complements?beneficiaryId=620903.6&donorId=479437.5&typeCategory=patricGenomeIds

    """
    beneficiaryId = request.args.get("beneficiaryId")
    donorId       = request.args.get("donorId")
    typeCategory  = request.args.get("typeCategory")

    try:
        compls = seed_base_response(
            beneficiaryId, donorId, typeCategory, all_seed_compls
        )
    except Exception as e:
        return jsonify({"error": str(e)}), 500

    return val_compl_res(
        compls, beneficiaryId, donorId, typeCategory, "seed complementarities"
    )


@app.route("/seed-scores", methods=["GET"])
def seed_scores():
    """
    For a pair of NCBI Taxonomy or GTDB representative genomes IDs or PATRIC IDs
    returns the seed score for all combinations between
    the beneficiary and the donor PATRIC ids found on microbetagDB.

    Examples:

        http://localhost:1337/seed-scores?beneficiaryId=86180&donorId=2184738&typeCategory=ncbiTaxonomyIds
        http://localhost:1337/seed-scores?beneficiaryId=86180.3&donorId=2184738.11&typeCategory=patricGenomeIds
        http://localhost:1337/seed-scores?beneficiaryId=GCF_002980625.1&donorId=GCF_011366425.1&typeCategory=ncbiGenomeIds

    """
    beneficiaryId = request.args.get("beneficiaryId")
    donorId       = request.args.get("donorId")
    typeCategory  = request.args.get("typeCategory")

    scores = seed_base_response(
        beneficiaryId, donorId, typeCategory, all_seed_scores
    )

    return val_compl_res(
        scores, beneficiaryId, donorId, typeCategory, "seed scores"
    )


# -----------------------------------------------------------------------------
# ------------------------          ON-THE-FLY        -------------------------
# -----------------------------------------------------------------------------

@app.route("/upload-abundance-table-dev", methods=["GET", "POST"])
@with_traceback_on_errors()  # noqa: F405
def upload_dev():

    if request.method == "POST":

        # Generate a unique identifier for the request based on its content; an integer
        request_id = hash(frozenset(request.form.items()))

        # Check if the request has already been processed recently (within 20 seconds; set in the cache above)
        if cache.get(request_id):

            # Return the cached response for the duplicate request
            response = cache.get(request_id)
            return jsonify({"status": "success", "data": response})

        # Init vars
        user_arg_vars = UserArgsVars()    # noqa: F405
        mtg_app_vars  = MicrobetagVars()  # noqa: F405

        # Parse user's input
        users_input = request.get_json()
        (
            args,
            out_dir,
            abd_data_df

        ) = args_and_abd_data(users_input)  # noqa: F405

        # Check whether the analysis is possible through the server
        check1000(abd_data_df, users_input)  # noqa: F405

        # Get user's input files in the required format
        (
            args["abd_table_file"],
            args["edge_list_file"],
            args["metadata_file"]

        ) = build_data_files(users_input, out_dir, abd_data_df)  # noqa: F405

        # FlashWeave (no network given), FAPROTAX and manta: the network workflow of the pipeline
        pipeline_outputs = None
        if pipeline_available() and (
            not args["edge_list_file"] or args[user_arg_vars.FAPROTAX] or args[user_arg_vars.NET_CLUSTER]
        ):
            try:
                pipeline_outputs = run_network_workflow(
                    out_dir,
                    abundance_table = args["abd_table_file"],
                    network         = args["edge_list_file"],
                    metadata        = args["metadata_file"],
                    faprotax        = bool(args[user_arg_vars.FAPROTAX]),
                    net_cluster     = bool(args[user_arg_vars.NET_CLUSTER]),
                    flashweave_args = {
                        "heterogeneous": args[user_arg_vars.HETEROGENEOUS],
                        "sensitive"    : args[user_arg_vars.SENSITIVE],
                        "max_k"        : args[user_arg_vars.MAX_K],
                    },
                )
            except Exception:
                logger.error(traceback.format_exc())
                raise RuntimeError(
                    "microbetag could not run FlashWeave, FAPROTAX or manta on your data. "
                    "Please, check your abundance table and network."
                )

        # Load config.yml based on user's args
        try:
            conf_modified_yaml = sync_yaml(out_dir, args, pipeline_outputs)  # noqa: F405

        except Exception:
            raise ArgumentMissing(
                "Could not build the updated YAML configuration file."
            )  # noqa: F405

        # Get microbeag config
        with open(conf_modified_yaml, "r") as yaml_file:
            init_yaml = yaml.safe_load(yaml_file)

        try:
            config = Config(init_yaml, conf_modified_yaml)
        except Exception:
            error_msg = "Configuration file failed to build an instance of the microbetag.Config class."
            raise FileFormatError(error_msg)  # noqa: F405

        config.db_config = login  # noqa: F405

        # Abundance data df
        data = prep_abundance_table(config)  # noqa: F405
        try:
            if lineage_matching():
                # A database of the user's own genomes: the closest genomes of each taxon, from its lineage
                config.otf_seq_tax_df, config.repr_genomes_present = map_seq_to_catalogue(
                    abd_tab   = data,
                    tax_col   = mtg_app_vars.TAX_COL,
                    seqId     = mtg_app_vars.SEQ_COL,
                    tax_delim = args[user_arg_vars.DELIMITER],
                )
            else:
                config.otf_seq_tax_df, config.repr_genomes_present = map_seq_to_ncbi_tax_level_and_id(
                    abd_tab    = data,
                    tax_col    = mtg_app_vars.TAX_COL,
                    seqId      = mtg_app_vars.SEQ_COL,
                    tax_scheme = args[user_arg_vars.TAXONOMY],
                    tax_delim  = args[user_arg_vars.DELIMITER],
                    mappings   = Mappings()  # noqa: F405
                )

        except Exception:
            # Catch and return any error that occurs
            error_msg = (
                "microbetag could not map taxonomies provided to their corresponding NCBI Tax IDs"
                "and/or GTDB representative genomes."
                "Please, check on the format of your taxonomy scheme, the delimiter provided etc."
            )
            logger.error(error_msg + "\n" + traceback.format_exc())
            raise FileFormatError(error_msg)  # noqa: F405

        taxonomy_file = os.path.join(out_dir, "map_tax_output.csv")
        config.otf_seq_tax_df.to_csv(taxonomy_file)

        # NOTE (Haris Zafeiropoulos, 2025-05-16):
        # Pass Julia running instance to the configuration instance
        # This saves time for initiating a new one every time
        config.jl = jl

        # The net_cx2 is a ndex2.cx2.CX2Network object
        net_cx2 = run_microbetag(config)

        # Return as response (JSONArray in Java, i.e. list[dict, dict, ..])
        return net_cx2.to_cx2()

    else:
        raise RuntimeError("Only POST queries supported for the whole pipeline run.")


if __name__ == "__main__":
    app.run()
