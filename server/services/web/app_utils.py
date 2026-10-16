import os
import yaml
import json
import traceback
from functools import wraps
from datetime import datetime

import pandas as pd

from dotenv import load_dotenv
from flask import Flask, jsonify

from microbetag import (
    Config, mtg_logger, __version__,
    db
)

__all__ = [
    'Mappings', 'UserArgsVars', 'MicrobetagVars', 'ArgumentMissing', 'NotFoundError', 'FileFormatError',
    'init_app', 'build_data_files', 'sync_yaml', 'args_and_abd_data', 'check1000', 'prep_abundance_table',
    'with_traceback_on_errors',
]


_logger_ = mtg_logger(__file__)


# -----------------
# Global variables
# -----------------
app_dir = os.path.dirname(os.path.abspath(__file__))


# Database credential
def db_credentials():

    __script_dir__ = os.path.dirname(os.path.abspath(__file__))

    try:
        load_dotenv()
        db_config = {
            "user"    : os.environ["USER_NAME"],
            "password": os.environ["PASSWORD"],
            "host"    : os.environ["HOST"],
            "port"    : os.environ["PORT"],
            "database": os.environ["DB_NAME"],
        }

        return db_config

    except Exception:

        raise RuntimeError("Database credentials could not be retrieved.")


# Fixed values for variables used in files microbetag builds itself
class MicrobetagVars:
    VERSION         = __version__.replace(".","")
    UPLOAD_FOLDER   = "/tmp"
    ABD_TABLE_DELIM = "\t"
    TAX_COL         = "taxonomy"
    TAX_DELIM       = ";"
    SEQ_COL         = "seqId"
    COM_CHAR        = "#"


# Parameters to be provided by the user, the name of which needs to be synced with MGG
# In case something is missing, gets default values
class UserArgsVars:
    """
    This class is being used to get attribute names, not actual values, across the project.
    The validate feature allows user's input to be updated with default values if something is missing.

    The names we give to the attributes here, need to be exactly the same as provided by the MGG.
    See https://github.com/hariszaf/MGG/blob/MGG/src/main/java/be/kuleuven/mgG/internal/tasks/SendDataToServerTask.java
    """

    PARAMETERS     = "parameters"

    INPUT_TYPE     = "input_type"
    ABUNDANCE_DATA = "abundance_data"
    TAXONOMY       = "taxonomy"
    DELIMITER      = "delimiter"

    NETWORK        = "network"
    METADATA       = "metadata"
    PHEN_TRAITS    = "phen_traits"
    FAPROTAX       = "faprotax"
    PATH_COMPL     = "pathway_complementarity"
    SEED_COMPL     = "seed_complementarity"
    SENSITIVE      = "sensitive"
    MAX_K          = "max_k"
    HETEROGENEOUS  = "heterogeneous"
    NET_CLUSTER    = "network_clustering"
    GET_CHILDREN   = "get_children"

    # Required top-level fields with optional allowed values
    _required = {
        INPUT_TYPE: ["abundance table", "network"],
        ABUNDANCE_DATA: {},  # Structure validated separately
        DELIMITER: [";", "|", "__", "_"],
        TAXONOMY: ["GTDB", "Silva", "microbetag_prep", "other"],
    }  # NOTE (Haris Zafeiropoulos, 2025-05-13): PARAMETERS is required but validated separately

    _extra_data = {
        NETWORK: None,
        METADATA: None,
    }
    _extra_params = {
        PHEN_TRAITS  : False,
        FAPROTAX     : True,
        PATH_COMPL   : False,
        SEED_COMPL   : False,
        SENSITIVE    : True,
        MAX_K        : 3,
        HETEROGENEOUS: False,
        NET_CLUSTER  : False,
        GET_CHILDREN : False
    }

    @classmethod
    def get_optional_fields(cls):
        return [
            getattr(cls, attr)
            for attr in dir(cls)
            if
            not attr.startswith("_") and isinstance(getattr(cls, attr), str) and getattr(cls, attr) not in cls._required
        ]

    @classmethod
    def validate(cls, input_data: dict):
        """
        input_data is the dictionary with the user's data and arguments
        This method aims to add entries in the dictionary that may be missing and provide default values.
        """

        if cls.PARAMETERS not in input_data:
            raise ArgumentMissing(
                f"Missing required argument '{cls.PARAMETERS}'. "
                "Please provide the 'parameters' dictionary with steps to go for."
            )

        args = input_data[user_vars.PARAMETERS]

        input_data[user_vars.PARAMETERS] = cls.split(args)

        for key, allowed_values in cls._required.items():

            in_top    = key in input_data
            in_params = key in input_data.get(cls.PARAMETERS, {})

            if not (in_top or in_params):

                if isinstance(allowed_values, list) and allowed_values:
                    valid_values = ", ".join(allowed_values)
                    hint         = f"{key} can take values: {valid_values}"

                else:
                    hint = "This field is required."

                raise ArgumentMissing(f"Missing required argument '{key}'. {hint}")

        for key, default_value in cls._extra_data.items():
            if key not in input_data:
                input_data[key] = default_value

        for key, default_value in cls._extra_params.items():
            if key not in input_data[cls.PARAMETERS]:
                input_data[cls.PARAMETERS][key] = default_value

        return input_data

    @classmethod
    def split(cls, orig_args: dict | list):

        if isinstance(orig_args, list):

            # If user does not set an argument MGG will set it as false
            args = {}
            for arg in orig_args:

                argument, value = arg.split(":")
                if value == "true":
                    evalue = True
                elif value == "false":
                    evalue = False
                else:
                    evalue = value

                args[argument] = evalue

        elif isinstance(orig_args, dict):
            args = orig_args

        else:
            raise ArgumentMissing(
                "Arguments need to be provided as a dictionary, as value to `parameters` key."
            )

        return args


# Paths to mapping files
class Mappings:
    """
    Paths to mapping files to be used
    """
    TAXONOMY          = os.path.join(app_dir, "mappings/taxonomy")
    SPECIES_NCBI_IDS  = os.path.join(TAXONOMY, "species2ncbiId.tsv")
    GENERA_NCBI_IDS   = os.path.join(TAXONOMY, "genera2ncbiId.tsv")
    FAMILIES_NCBI_IDS = os.path.join(TAXONOMY, "families2ncbiId.tsv")

    # The seed sets of the genomes are in the database (seedSets, nonSeedSets): read per network
    SEED_KO_MO          = os.path.join(app_dir, "mappings/kegg_seeds/seedId_keggId_module.tsv")


mtg_vars  = MicrobetagVars()
user_vars = UserArgsVars()
mappings  = Mappings()

# NOTE (Haris Zafeiropoulos, 2025-06-08): is this necessary ?
db.DB_CREDENTIALS = db_credentials()


# Errors - cases and formatting to be returned on the user-side
class NotFoundError(Exception):
    """Custom exception for genome lookup failures."""
    pass


class ArgumentMissing(Exception):
    """Custom excepion for missing argument"""
    pass


class FileFormatError(Exception):
    pass


def with_traceback_on_errors(*exception_classes):
    """
    Decorator that catches specified exceptions raised in a function and re-raises them
    with the full traceback included in the error message.

    If no exception classes are provided, it defaults to catching all exceptions.
    It is useful to avoid repeating the try-except pattern everywhere.

    Arguments:
        *exception_classes : Exception
            One or more exception classes to catch (e.g., ValueError, FileFormatError).
            If omitted, the decorator catches all exceptions.

    Returns:
        callable
            A wrapped function that raises exceptions with appended traceback information.

    Example:
        @with_traceback_on_errors(FileFormatError, ValueError)
        def risky_function():
            raise ValueError("Something went wrong")
    """
    # Default to all exceptions if none provided
    if not exception_classes:
        exception_classes = (Exception,)

    def decorator(func):

        @wraps(func)
        def wrapper(*args, **kwargs):

            try:
                return func(*args, **kwargs)
            except exception_classes as e:
                tb = traceback.format_exc()
                return jsonify({"error": str(e), "traceback": tb}), 500

        return wrapper

    return decorator

# ---------
# Init app
# ---------

def init_app():
    """
    Function to fire the app; it is performed only once.
    """
    app                              = Flask(__name__)
    app.debug                        = True
    app.config["UPLOAD_FOLDER"]      = mtg_vars.UPLOAD_FOLDER
    app.config["MAX_CONTENT_LENGTH"] = 16 * 1000 * 1000

    load_dotenv()
    app.secret_key = os.environ["FLASK_SECRET_KEY"]

    # FlashWeave in the server's Julia (up to 1.0.4); from 2.0, the network is an input
    jl = None
    if os.environ.get("FLASHWEAVE", "true").lower() == "true":
        from julia.api import Julia

        jl = Julia(compiled_modules=False)
        jl.using("FlashWeave")

    db = db_credentials()  # noqa: F405

    return app, db, jl


# --------------
# For the on-the-fly pipeline
# --------------

def args_and_abd_data(input_data):
    """
    Parse user's arguments to get a data frame with their abundance data 
    and a dictionary with their arguments/parameters on the MGG Parameters settings.
    """
    # Check whether all required user arguments are there
    # updated_data = user_vars.validate(input_data)
    user_vars.validate(input_data)

    # Create output dir for the session
    out_dir = _onthefly_outdir()

    try:
        # Empty lines (e.g. a blank last line of the user's table) are not taxa
        rows     = [row for row in input_data[user_vars.ABUNDANCE_DATA] if row and any(str(v).strip() for v in row)]
        abd_data = pd.DataFrame(rows)

    except Exception:
        error_msg = (
            f"No '{user_vars.ABUNDANCE_DATA}' provided."
            f"{user_vars.ABUNDANCE_DATA} is required and it is supposed to be provided as a list of lists, "
            "where each entry stands for a line and its entries are the values of your abundance table. "
            "For example:\n"
            "['ASV1076', '0', '0', '0', 'my_cool_speciies']"
        )
        raise ArgumentMissing(error_msg)

    args = input_data[user_vars.PARAMETERS]

    return args, out_dir, abd_data


def check1000(data: pd.DataFrame, input_data: dict) -> None | FileFormatError:
    """
    Checks whether abundance data to be analysed have up to 1000 sequence ids
    or, in case they have more, a network is also provided and the taxonomy scheme is not 'Other'.

    Args:
        data: A data frame with the abundance data to be analysed
        input_data: A dictionary with the user's input from the query
    """
    if data.shape[0] > 1000:

        accepted_tax_schemes = ["GTDB", "Silva", "microbetag_prep"]

        error_msg = ""

        if (user_vars.NETWORK not in input_data):

            error_msg += (
                "Your abundance table consists of more than 1000 sequences. "
                "Please infer a co-occurrence network first and try again microbetag providing the network too."
                "The microbetag_prep acompanying tool can be of help for network inference with FlashWeave."
                "https://github.com/hariszaf/microbetag/tree/preprocess"
            )

        if input_data[user_vars.PARAMETERS][user_vars.TAXONOMY] not in accepted_tax_schemes:
            valid_taxonomies = ", ".join(x for x in accepted_tax_schemes)
            error_msg += (
                "When running microbetag with an abundance table of more than 1000 sequences,"
                f"you need to make sure your taxonomy scheme is one of the following: {valid_taxonomies}"
            )

        if input_data[user_vars.PARAMETERS][user_vars.NET_CLUSTER]:
            error_msg += (
                "On the fly network clustering with manta is only available for abundance tables up to 1000 taxa."
                "Consider running microbetag without the network clustering option and then run manta individually."
                "Then you can add columns in the nodes table on Cytoscape with the cluster assigned to each node."
                "You may check here for more: https://microbetag.readthedocs.io/en/latest/advanced_use/large_data.html"
            )

        if error_msg != "":
            _logger_.error(error_msg)
            raise FileFormatError(error_msg)

    return None


def build_data_files(input_data: dict, out_dir: str, data: pd.DataFrame) -> tuple[str, str, str | None]:
    """
    Build data files to be able to be used in the app.

    Returns:
        A tuple including:
            - abundance_table_file: Filepath to the .csv abundance table on microbetag server
            - edge_list_file: Filepath to edgelist
            - metadata_file: Filepath to metadata file (optional)
    """
    if input_data[user_vars.PARAMETERS][user_vars.INPUT_TYPE] == "network" and input_data[user_vars.NETWORK]:

        try:
            # Skip header line
            # NOTE (Haris Zafeiropoulos, 2025-05-14):
            # The first row will be the headers, we keep only the list elements with the actual associatittons
            network      = pd.DataFrame(input_data[user_vars.NETWORK][1:])
            network_path = os.path.join(out_dir, "edgelist.tsv")

            network.to_csv(network_path, sep="\t", header=False, index=False)

            edge_list_file = network_path

        except Exception:
            error_msg = (
                "Network file was not provided properly and could not be exported as .csv in microbetag server."
            )
            _logger_.error(error_msg)
            raise FileFormatError(error_msg)
    else:
        edge_list_file = None

    if input_data[user_vars.METADATA]:

        try:
            metadata      = pd.DataFrame(input_data[user_vars.METADATA])
            metadata_path = os.path.join(out_dir, "metadata.tsv")
            metadata.to_csv(metadata_path, sep="\t", header=False, index=False)

            metadata_file = metadata_path

        except Exception:
            error_msg = "Unable to load/save metadata file on the server side."
            raise FileFormatError(error_msg)
    else:
        metadata_file = None

    # Write the user's abundance table as a .tsv file
    abundance_table_file = os.path.join(out_dir, "abundance_table.tsv")
    data.iloc[0, 0]      = mtg_vars.SEQ_COL
    data.iloc[0, -1]     = mtg_vars.TAX_COL
    data.to_csv(abundance_table_file, sep="\t", header=False, index=False)

    return abundance_table_file, edge_list_file, metadata_file


def sync_yaml(out_dir: str, args: dict, pipeline_outputs: dict = None) -> str:
    """
    This function edits the base yaml file based on the user's arguments.
    It will save a new run-specific .yaml file in the output folder of the run, 
    which will be then used with the Config() class of microbetag to fire the pipeline.

    In MGG mgG/internal/tasks/SendDataToServerTask.java:
        parameters.add("INPUT:" + INPUT.getSelectedValue());

    Returns:
        Filepath to updated YAML configuration file
    """

    yaml_version = f"onthefly_v{mtg_vars.VERSION}.yml"

    with open(yaml_version) as f:
        base_yaml = yaml.safe_load(f)

    base_yaml["output_directory"]["dir_path"] = out_dir

    # User's data files as processed by build_data_files()
    base_yaml["abundance_table_file"]["file_path"] = os.path.join(out_dir, "abundance_table.tsv")
    base_yaml["edge_list"]["file_path"]            = args["edge_list_file"]
    base_yaml["metadata_file"]["file_path"]        = args["metadata_file"]

    # Taxonomy and flashweave arguments
    base_yaml["taxonomy"] = args[user_vars.TAXONOMY]

    if "flashweave_args" in base_yaml:
        base_yaml["flashweave_args"]["heterogeneous"]["value"] = args[user_vars.HETEROGENEOUS]
        base_yaml["flashweave_args"]["sensitive"]["value"]     = args[user_vars.SENSITIVE]
        base_yaml["flashweave_args"]["max_k"]["value"]         = args[user_vars.MAX_K]

    # Steps; from 2.0, FAPROTAX and manta run in the pipeline, so without it they stay off
    from app_pipeline import pipeline_available

    unavailable = []
    if "flashweave_args" not in base_yaml and not pipeline_available():
        unavailable = ["faprotax_annotation", "network_clustering"]
    for step, arg in [
        ("faprotax_annotation", user_vars.FAPROTAX),
        ("phenotrex_traits", user_vars.PHEN_TRAITS),
        ("pathway_complementarity", user_vars.PATH_COMPL),
        ("seed_complementarity", user_vars.SEED_COMPL),
        ("network_clustering", user_vars.NET_CLUSTER),
    ]:
        if step in unavailable:
            if args[arg]:
                _logger_.warning(f"{step} is not available on this server; skipped.")
            continue
        base_yaml.setdefault(step, {})["value"] = args[arg]

    # Outputs of the network workflow of the pipeline (app_pipeline.run_network_workflow)
    pipeline_outputs = pipeline_outputs or {}
    if "edge_list" in pipeline_outputs:
        base_yaml["edge_list"]["file_path"] = pipeline_outputs["edge_list"]
    if "faprotax_sub_tables" in pipeline_outputs:
        base_yaml["faprotax_sub_tables"] = {"dir_path": pipeline_outputs["faprotax_sub_tables"]}
    if "manta_net" in pipeline_outputs:
        base_yaml["prev_clustered_network"] = {"file_path": pipeline_outputs["manta_net"]}

    modified_yaml = os.path.join(out_dir, "config.yml")

    with open(modified_yaml, "w") as f:
        yaml.safe_dump(base_yaml, f, sort_keys=False)

    return modified_yaml


def prep_abundance_table(conf: Config):
    """
    Wrap function for building abundance table in microbetag-orientd format.

    Arguments:
        conf: Instance of the microbetag `Config` class, with the extra flags of the sync_yaml().

    Returns:
        abd_df (pd.DataFrame): A data frame with the user's bundance data the way FlashWeave, in 
            the microbetag-specific implementation. 
    """
    try:
        abundance_table = _read_abundance_table(conf.abundance_table)
    except Exception:
        error_msg = "Failed to read initial abundance data file during the preprocess step."
        raise FileFormatError(error_msg)

    if os.path.exists(conf.network):
        abundance_table["microbetag_id"] = abundance_table[mtg_vars.SEQ_COL]
        return abundance_table

    try:
        abd_df = _prepare_flashweave_format(abundance_table, conf.flashweave_abd_table)
    except Exception:
        error_msg = "Failed to prepare abundance table in the format required for FlashWeave to run."
        _logger_.error(error_msg)
        return FileFormatError(error_msg)  # jsonify({"status": "error", "message": error_msg}), 400

    return abd_df


def _onthefly_outdir() -> str:
    """Create output dir for the session: its timestamp, and a random suffix for requests of the same second"""
    import tempfile

    dt_string = datetime.now().strftime("%d-%m-%Y--%H.%M.%S")
    return tempfile.mkdtemp(prefix=f"{dt_string}-", dir=mtg_vars.UPLOAD_FOLDER)


def _read_abundance_table(filepath) -> pd.DataFrame | FileFormatError:
    n_skip = _count_comment_lines(filepath)

    try:

        df = pd.read_csv(filepath, sep="\t", skiprows=n_skip)
        if df.shape[1] < 2:
            raise FileFormatError("No valid records found.")
        return df

    except Exception as e:
        _logger_.error("Invalid OTU table format. Use tab- or comma-separated values.")
        raise FileFormatError(str(e)) from e


def _count_comment_lines(filepath) -> int:
    with open(filepath, "r") as f:
        for index, line in enumerate(f):
            if not line.startswith(mtg_vars.COM_CHAR):
                return index
    return 0


def _prepare_flashweave_format(df, out_file=None):
    """
    FlashWeave needs an abundance table where taxonomy column is skipped, numeric columns are integer,
    and not float, and, for microbetag, rows need to be taxa and columns the samples.
    """
    try:

        fw_table = df.drop(columns=mtg_vars.TAX_COL).copy()

        fw_table[
            fw_table.select_dtypes("float64").columns
        ] = fw_table.select_dtypes("float64").astype("int64")

        fw_table[mtg_vars.SEQ_COL] = fw_table[mtg_vars.SEQ_COL].astype(str)

        # Export df to .csv file
        fw_table.to_csv(out_file, sep="\t", index=False)

        df["microbetag_id"] = fw_table[mtg_vars.SEQ_COL]

    except Exception:

        error_msg = "Failed to prepare abundance table in the format required for FlashWeave to run."
        _logger_.error(error_msg)
        raise FileFormatError(error_msg)  # jsonify({"status": "error", "message": error_msg}), 400

    return df
