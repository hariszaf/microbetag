# microbetag : a software suite to annotate microbial co-occurrence networks

# Copyright (c) 2025 Haris Zafeiropoulos

# Licensed under GNU LGPL.3, see LICENCE file

"""
Utility functions to be used across the `microbetag` library.
"""

import os
import sys
import csv
import ast
import json
import copy
import shutil
import logging
import colorlog
import pandas as pd
from pathlib import Path
from multiprocessing import Pool, cpu_count
from typing import Any, Dict, List, Tuple, Set, Union

# Handling data related


def resolve_relative_path(base_dir: str, file_path: str) -> str:
    """
    Resolves a relative file path into an absolute file path based on a given base directory.

    This function processes a relative `file_path` (which may contain one or more `../` 
    segments) and resolves it into an absolute path by moving back the corresponding 
    number of directory levels from `base_dir`. It returns the resulting absolute file path.

    Parameters
    ----------
    base_dir : str
        The base directory from which to resolve the relative `file_path`. This should
        be an absolute path to a directory.

    file_path : str
        The relative file path to be resolved. It may contain `../` to navigate up the 
        directory hierarchy.

    Returns
    -------
    str
        The resolved absolute file path.

    Examples
    --------
    >>> resolve_relative_path("/home/user/docs", "../files/report.txt")
    '/home/user/files/report.txt'

    """

    return str(Path(base_dir).resolve().joinpath(file_path).resolve())


def resolve_file_path(base_dir: str, file_path: str) -> str:
    """
    Resolves a file path relative to a given base directory and returns the absolute file path.

    If the provided `file_path` is relative, it is resolved using the `base_dir`. The function
    handles absolute paths, user directory expansion (e.g., `~`), and relative paths (e.g., `../`).

    Parameters
    ----------
    base_dir : str
        The base directory to resolve relative paths from.

    file_path : str
        The file path to resolve. It can be absolute, relative, or use `~` for the home directory.

    Returns
    -------
    str
        The resolved absolute file path.

    Raises
    ------
    FileNotFoundError
        If the resolved file path does not exist.

    Examples
    --------
    >>> resolve_file_path("/home/user/docs", "~/file.txt")
    '/home/user/file.txt'
    """

    _logger_.info(f">> file_path: {file_path}")

    if file_path is None:
        return None  # Return None if the file path is None

    # Handle absolute paths and ~ expansion
    if file_path.startswith("/"):
        path = Path(file_path).resolve()

    elif file_path.startswith("~"):
        path = Path(os.path.expanduser(file_path)).resolve()

    elif file_path.startswith("../"):
        # Resolve relative path using the base_dir
        path = Path(resolve_relative_path(base_dir, file_path)).resolve()

    else:
        # Relative path with respect to base_dir
        path = (Path(base_dir) / file_path).resolve()

    # Check if the file exists
    if path.exists():
        return str(path)
    else:
        raise FileNotFoundError(f"File not found: {path}")


class SetEncoder(json.JSONEncoder):
    """
    Custom JSON encoder that handles serialization of Python sets.

    This encoder extends the functionality of the standard JSONEncoder to support
    serializing Python sets. JSON does not have a native representation for sets,
    so this encoder converts sets to lists before serializing them.

    Usage:
        When serializing data to JSON using json.dump() or json.dumps(), specify
        cls=SetEncoder to use this custom encoder.

    References:
        - json.JSONEncoder: https://docs.python.org/3/library/json.html#json.JSONEncoder
    """

    def default(self, obj):
        """
        Override the default method of JSONEncoder to handle serialization of sets.
        Notes:
            If the object is a set, it is converted to a list before serialization.
            Otherwise, the default behavior of JSONEncoder.default() is used.        Notes:
            If the object is a set, it is converted to a list before serialization.
            Otherwise, the default behavior of JSONEncoder.default() is used.
        """
        if isinstance(obj, set):
            return list(obj)
        return json.JSONEncoder.default(self, obj)


def mtg_logger(filename: str) -> logging.getLogger:
    """
    Creates and returns a configured logger instance. This logger:

    - Logs messages to stdout with colored formatting using `colorlog`
    - Avoids adding duplicate handlers if called multiple times
    - Uses the given `filename` as the logger's name
    - Logs messages with level INFO and above

    Arguments:
        script: The filename of the script where the logger will be applied to.

    Returns:
        The logger instance.
    """
    logger = logging.getLogger(filename)
    logger.setLevel(logging.INFO)
    if not logger.handlers:
        sh = logging.StreamHandler(sys.stdout)
        sh.setLevel(logging.INFO)
        formatter = colorlog.ColoredFormatter(
            "%(log_color)s%(asctime)s - %(name)s - %(levelname)s - %(message)s",
            log_colors={
                "DEBUG": "cyan",
                "INFO": "blue",
                "WARNING": "yellow",
                "ERROR": "red",
                "CRITICAL": "bold_red",
            },
        )
        sh.setFormatter(formatter)
        logger.addHandler(sh)
        logger.propagate = False  # Prevent duplicate logging
    return logger


def safe_literal_eval(value: Any):
    """
    Safely evaluates a string that may represent a Python literal (e.g., list, dict, int).

    This function attempts to parse a string using :class:`ast.literal_eval`, which only evaluates
    Python literals (e.g., strings, numbers, tuples, lists, dicts, booleans, and None),
    avoiding the security risks of `eval()`*.
    If `value` is not a string or if evaluation fails,
    the original value is returned unchanged.

    Args
    ----------
    value : Any
        The input to be evaluated. If it's a string that looks like a literal (e.g., "[1, 2]"), 
        it will be parsed. Otherwise, it's returned as is.

    Returns
    ----------
        The evaluated literal if successful, or the original value if evaluation fails.

    Examples
    ----------
    >>> safe_literal_eval("[1, 2, 3]")
    [1, 2, 3]

    >>> safe_literal_eval("{'a': 1}")
    {'a': 1}

    Note:
        * Security risks of `eval`:
        https://www.adventuresinmachinelearning.com/safe-and-secure-eval-in-python-how-to-minimize-security-risks/
    """
    try:
        # Attempt to evaluate the value if it's a string that looks like a list
        return ast.literal_eval(value) if isinstance(value, str) else value
    except (ValueError, SyntaxError):
        # If it's not a valid list string, return the original value
        return value


def flatten(list_of_lists: List) -> List:
    """
    Recursively flattens a nested list into a single-level list.

    This function handles arbitrarily nested lists and returns a new list
    containing all the leaf elements in the original order.

    Arguments:
        lst : A list that may contain other nested lists.

    Returns:
    -------
        A flat list containing all non-list elements in the original order.

    Examples:
    --------
    >>> flatten([1, [2, [3, 4]], 5])
    [1, 2, 3, 4, 5]
    """
    if len(list_of_lists) == 0:
        return list_of_lists
    if isinstance(list_of_lists[0], list):
        return flatten(list_of_lists[0]) + flatten(list_of_lists[1:])
    return list_of_lists[:1] + flatten(list_of_lists[1:])


def flatten_list(lst: List, flat_list: List = None) -> set:
    """
    Recursively flattens a nested list and returns a set of unique elements.

    This function traverses all nested lists and collects elements into a set,
    removing any duplicates. The final result is unordered.

    Parameters:
    ----------
    lst: A list that may contain other nested lists.
    flat_list : Optional list. Used internally during recursion. Should not be set manually.

    Returns:
    -------
        A set containing all unique elements from the nested list.

    Examples:
    --------
    >>> flatten_list([1, [2, [2, 3]], 4, 1])
    {1, 2, 3, 4}
    """
    if flat_list is None:
        flat_list = []

    for i in lst:
        if isinstance(i, list):
            flatten_list(i, flat_list)
        else:
            flat_list.append(i)

    return set(flat_list)


def load_merged_ko_file(merged_ko: str) -> pd.DataFrame:
    """
    Load the 3-columns KEGG annotations file as built from the merge_ko()

    Input:
        merged_ko: path to 3-columns output file of the merge_ko()

    Returns:
        pivot_df: a presence-absence (1/0) df where KOs are the rows and bin_ids the columns
    """
    p = Path(merged_ko)

    # Handle gzipped file
    if p.suffix == ".gz":
        decompressed = p.with_suffix("")  # remove .gz
        if not decompressed.exists():
            import gzip
            _logger_.info(f"Decompressing {p.name} -> {decompressed.name}")
            with gzip.open(p, 'rb') as f_in, open(decompressed, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        merged_ko = decompressed
    else:
        # Compress if uncompressed version found and no .gz exists
        gz_path = p.with_suffix(p.suffix + ".gz")
        if not gz_path.exists():
            # Keep the original file: it is read right below
            import gzip
            with open(p, "rb") as f_in, gzip.open(gz_path, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        merged_ko = p

    # Load the (decompressed) file
    df = pd.read_csv(merged_ko, sep="\t")

    # Build presence/absence pivot table
    bin_id, _, ko = df.columns[:3]  # The second column of the file corresponds to the contig id.
    unique   = df.drop_duplicates().assign(presence=1)
    pivot_df = unique.pivot_table(index=ko, columns=bin_id, values="presence", fill_value=0)

    return pivot_df  # keep one | used to alse return the bins_kos


def convert_to_json_serializable(obj: Any) -> Any:
    """
    Recursively serializes entries of an object
    A set is converted to a list, a list is flattened to its items
    and a dictionary keeps its key and their values get serialized.

    Note:
        This is essential step both for allowing a jsonified response and to be able
        to dump a dictionary as a JSON file.
    """
    if isinstance(obj, (int, float, str, bool, type(None))):
        return obj

    elif isinstance(obj, set):
        return list(obj)

    elif isinstance(obj, list):
        return [convert_to_json_serializable(item) for item in obj]

    elif isinstance(obj, dict):
        new_dict = {}
        for key, value in obj.items():
            if not isinstance(key, (str, int, float, bool, type(None))):
                key = str(key)  # or use "|".join(key) if you want to preserve tuple structure better
            new_dict[key] = convert_to_json_serializable(value)
        return new_dict

    else:
        try:
            return json.dumps(obj)
        except TypeError:
            return str(obj)


def _extend_for_pair(args):
    beneficiary_bin, potential_donor, compls, descrps, pc_percent = args
    local_out = {}

    if not compls:
        return (beneficiary_bin, potential_donor, local_out)

    for compl in compls:
        module_id   = compl[0][3:] if compl[0].startswith("md") else compl[0]
        kos_to_get  = compl[1]
        complet_alt = compl[2]

        if len(kos_to_get) / len(complet_alt) > pc_percent:
            continue

        compl_str = [x if isinstance(x, str) else ";".join(x) for x in compl[1:]]
        triplet   = descrps[descrps["moduleId"] == module_id].values.tolist()[0]

        local_out[len(local_out)] = triplet + compl_str

    return (beneficiary_bin, potential_donor, local_out)


def extend_complements(
    complements_json: str,
    descrps_path: str,
    pc_percent: int,
    pc_dir: str,
    n_workers: int = None
):
    """
    Extends pathway complement annotations based on given settings and descriptions.

    Parameters:
        - complements_json: Path to the complements JSON file.
        - descrps_path: Path to the KEGG MODULES description file.
        - pc_percent: Maximum allowable percentage of required KOs that must be present.
        - pc_dir: Directory to save the extended complements JSON file.
        complements_dict (dict): Dictionary of complements loaded from a JSON file.
        descrps_path (str): Path to the module descriptions file (tab-separated file with no header).

    Returns:
        A dictionary with pathway complementarities to be assigned in the MGG format

    Note:
        Here we build the `pathway_complements_extended.json` a JSON file with the dictionary returned
    """

    if n_workers is None:
        n_workers = max(cpu_count() - 1, 1)

    descrps         = pd.read_csv(descrps_path, sep="\t", header=None)
    descrps.columns = ["category", "moduleId", "description"]
    descrps         = descrps[["moduleId", "description", "category"]]

    out_json = os.path.join(pc_dir, "pathway_complements_extended.json")
    if os.path.exists(out_json):
        return json.load(open(out_json))

    complements_dict = json.load(open(complements_json))
    complements_dict_ext = copy.deepcopy(complements_dict)

    tasks = []
    for beneficiary_bin, donors in complements_dict.items():
        for donor, compls in donors.items():
            tasks.append((beneficiary_bin, donor, compls, descrps, pc_percent))

    with Pool(n_workers) as pool:
        results = pool.map(_extend_for_pair, tasks)

    # merge
    for bin, donor, data in results:
        complements_dict_ext[bin][donor] = data

    json.dump(complements_dict_ext, open(out_json, "w"))
    return complements_dict_ext


def extend_faprotax(faprotax_sub_tables, sequence_id_column_name) -> Tuple[dict[str, list], list[str]]:
    """
    Parses the sub tables of the faprotax analysis
    to assign the biological processes related to each sequence id

    Returns:
        A tuple consisting of:
        - bin_faprotax_traits: A dictionary with the sequence id as key and a list of FAPROTAX trais a value
        - faprotax_traits: A list with the unique set of the FAPROTAX traits found across all taxa of the study
    """
    bin_faprotax_traits = {}

    fapro_sub_tables = [
        os.path.join(faprotax_sub_tables, file)
        for file in os.listdir(faprotax_sub_tables)
    ]

    for file in fapro_sub_tables:
        _logger_.info(f"File: {file}")
        # NOTE (Haris Zafeiropoulos, 2025-05-20):
        # We replace '_' with a space for user's convenience in the MGG
        # Also, this needs to be synced with the MGG.MUtils code for the grouping in the node panel
        trait_name, _ = os.path.splitext(os.path.basename(file))
        # FAPROTAX writes the sub-tables with the abundance table's delimiter; line 1 is a comment
        with open(file) as f:
            f.readline()
            sep = "\t" if "\t" in f.readline() else ","
        trait         = pd.read_csv(file, sep=sep, skiprows=1)

        bins_with_trait = trait[sequence_id_column_name].dropna()

        for bin_id in bins_with_trait:
            bin_faprotax_traits.setdefault(bin_id, []).append(trait_name.replace("_", " "))

    faprotax_traits = list(flatten_list(bin_faprotax_traits.values()))

    return bin_faprotax_traits, faprotax_traits


def load_phenotypic_traits(phen_outdir) -> Tuple[Dict[str, Dict[str, Union[str, float]]], Set[str]]:
    """
    Load phenotrex-based trait files and assignm them per genome.

    Returns:
        A tuple consisting of:
            - bin_phen_traits: A dictionary with genome id as key and a dictionary as value, 
                                with each phenotrex-based trait as value and their presence/absence
                                and scores as value
        phentraits: A set with the traits presentt

    Note:
        Example of a `bin_phen_traits`:
        ```
        bin_phen_traits[bin_id][trait_name] = {
            "presence": case["Trait present"],
            "confidence": case["Confidence"],
        }
        ```
    """
    bin_phen_traits = {}
    phentraits      = set()

    prediction_files = [
        os.path.join(phen_outdir, file)
        for file in os.listdir(phen_outdir)
    ]

    for file in prediction_files:

        if os.path.getsize(file) == 0:
            continue

        trait          = pd.read_csv(file, sep="\t", skiprows=1)
        trait_name     = os.path.basename(file).split(".prediction.tsv")[0]
        trait_filtered = trait[trait["Trait present"].notna()]
        trait_dict     = trait_filtered.to_dict(orient="records")

        for case in trait_dict:
            bin_id, _ = os.path.splitext(case["Identifier"])
            if bin_id not in bin_phen_traits:
                bin_phen_traits[bin_id] = {}
            phentraits.add(trait_name)
            bin_phen_traits[bin_id][trait_name] = {
                "presence": case["Trait present"],
                "confidence": case["Confidence"],
            }

    return bin_phen_traits, phentraits


from typing import Optional


def detect_separator(file_path: str) -> Optional[str]:
    """
    Detects the separator used in a text file.
    Tries csv.Sniffer first, falls back to common delimiters.

    Returns:
        Detected delimiter as a string or None if undetectable.
    """

    # Common delimiters to check if csv.Sniffer fails
    fallback_delimiters = ['\t', ',', ';', ' ']

    try:
        with open(file_path, "r", encoding="utf-8") as f:
            f.seek(0, 2)
            file_size = f.tell()

            if file_size == 0:
                raise ValueError("File is empty")

            # Use a simpler and safer sampling logic
            sample_size = max(min(file_size, int(file_size * 0.1)), 2000)
            f.seek(0)
            sample = f.read(sample_size)

            if not sample.strip():
                raise ValueError("File contains only whitespace")

            # Try Sniffer first
            try:
                sniffer = csv.Sniffer()
                dialect = sniffer.sniff(sample)
                return dialect.delimiter
            except Exception:
                pass  # Fall through to manual detection

            # Manual heuristic: test fallback delimiters
            first_line = sample.splitlines()[0]
            for delim in fallback_delimiters:
                if delim in first_line:
                    return delim

            # Last resort: try Sniffer again on smaller chunk
            try:
                dialect = csv.Sniffer().sniff(sample[:1000])
                return dialect.delimiter
            except Exception:
                return None

    except Exception as exc:
        raise TypeError(
            f"Could not detect delimiter for file: {file_path}. "
            f"Reason: {exc}"
        ) from exc


def find_three_column_format(file_path: str, delimiter: str) -> tuple[int, Union[None, int]]:
    """
    Checks if a file is in a three-column format and whether the third column is a float.
    If not, it skips row and goes to the next one checking for the 3-colummn format.
    Once met, it returns the line number, if that is neve the case raises an Exception.

    Args:
        file_path (str): Path to the file to be checked.
        delimiter (str): The delimiter used to separate columns (e.g., '\t' for tab-separated values).

    Returns:
        tuple: (line_number, None) if the third column is a float, (line_number, 0) otherwise.
    """
    with open(file_path, "r") as f:
        for line_num, line in enumerate(f, start=1):
            columns = line.strip().split(delimiter)
            if len(columns) == 3:
                try:
                    float(columns[-1])
                    return line_num, None
                except (ValueError, TypeError):
                    pass
    raise TypeError(
        f"The network file {file_path} is not in the 3-columns format required with a numeric weight column."
    )


_logger_ = mtg_logger(__name__)

