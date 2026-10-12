import argparse
import gzip
import json
import os
import re
import time

# SQLAlchemy
import urllib.parse

import pandas as pd
from sqlalchemy import create_engine, text

on_docker = "/var/lib/mysql-files"
MAPPINGS = (
    os.path.join(on_docker, "mappings")
    if os.path.exists(os.path.join(on_docker, "mappings"))
    else "data/mappings"
)
ANNOTATIONS = (
    os.path.join(on_docker, "annotations")
    if os.path.exists(os.path.join(on_docker, "annotations"))
    else "data/annotations"
)


class Credentials:
    def __init__(self, USER_NAME, PASSWORD, HOST, DB_NAME, PORT):
        self.USER_NAME = USER_NAME
        self.PASSWORD = PASSWORD
        self.HOST = HOST
        self.DB_NAME = DB_NAME
        self.PORT = PORT


# Create enging with mysql.connector
def create_cursor(credentials):
    """
    [DEPRECATED]
    Creates a connection to the microbetagDB; gets values from the variables.py.
    """
    import mysql.connector

    cnx = mysql.connector.connect(
        user=credentials.USER_NAME,
        password=credentials.PASSWORD,
        host=credentials.HOST,
        database=credentials.DB_NAME,
    )
    cnx.get_warnings = True
    cursor = cnx.cursor()
    return cnx, cursor


# Create an SQLAlchemy engine
# NOTE (Haris Zafeiropoulos, 2025-03-27): Needs to got with urllib.parse if password has special characters
def create_sqlalchemy_engine(credentials):
    PASSWORD_LIB = urllib.parse.quote(credentials.PASSWORD)
    engine = create_engine(
        f"mysql+pymysql://{credentials.USER_NAME}:{PASSWORD_LIB}@{credentials.HOST}:{credentials.PORT}/{credentials.DB_NAME}"
    )
    return engine


def load_db_credentials(database):
    with open(".env_dev.json", "r") as f:
        credentials = json.load(f)

    USER_NAME, PASSWORD, HOST, DB_NAME, PORT = credentials[database].values()

    return Credentials(USER_NAME, PASSWORD, HOST, DB_NAME, PORT)


def get_table_columns(engine, table_name):
    with engine.connect() as connection:
        raw_connection = connection.connection
        cursor = raw_connection.cursor()
        cursor.execute(f"SELECT * FROM {table_name} LIMIT 0")
        column_names = [desc[0] for desc in cursor.description]
    return column_names


# Load TSV to SQL
def load_tsv_to_sql(
    engine,
    file,
    table_name,
    encoding="utf-8",
    header=None,
    if_exists="append",
    trim_suffix=False,
    column=None,
    escape_char=None,
    first_column_index=False,
):
    """
    Loads a TSV file into a SQL database table using an SQLAlchemy engine.

    Parameters:
    engine (sqlalchemy.engine.base.Engine): The SQLAlchemy engine connected to the database.
    file (str): The path to the TSV file to be loaded. Supports both plain and gzip-compressed files.
    table_name (str): The name of the database table where the TSV data will be loaded.
    if_exists (str, optional): The behavior when the table already exists. Defaults to "append".
                               Options include 'fail', 'replace', or 'append'.

    """
    # Get column names
    column_names = get_table_columns(engine, table_name)

    if file.endswith(".gz"):
        with gzip.open(file, "rt", encoding=encoding) as f:
            df = pd.read_csv(f, sep="\t", header=header)
    else:
        df = pd.read_csv(file, sep="\t", header=header)

    if first_column_index:
        df = df.iloc[:, 1:]

    df = df.apply(
        lambda col: (
            col.str.replace(r"[\r\n]", "", regex=True) if col.dtype == "object" else col
        )
    )

    df.columns = column_names

    if trim_suffix:
        # df[column].str.replace(r"\.[^.]*$", "", regex=True)
        escaped_char = re.escape(escape_char)
        df[column] = df[column].str.replace(
            rf"{escaped_char}[^{escaped_char}]*$", "", regex=True
        )

    df.to_sql(table_name, con=engine, if_exists=if_exists, index=False)


# Seed tables
def load_seed_scores(engine, file):
    """
    Load seed scores to the seedScores table.
    """
    if file.endswith(".gz"):
        with gzip.open(file, "rt", encoding="utf-8") as f:
            df = pd.read_csv(f, sep="\t", header=None)
    else:
        df = pd.read_csv(file, sep="\t", header=None)

    print(f"File {file} was loaded as dataframe.")

    if "competitionScore" not in df.columns:
        df.columns = [
            "patricGenomeA",
            "patricGenomeB",
            "competitionScore",
            "complementaritScore",
        ]

    df.to_sql("seedScores", con=engine, if_exists="append", index=False)

    print(f"File {file} was successfully loaded to the seed scores table.")


def create_index_for_seed_scores(engine):
    with engine.connect() as connection:
        # Create an index on the 'patricId' and 'genomeId' columns
        create_index_query = """
            CREATE INDEX genome_pair ON seedScores(patricGenomeA, patricGenomeB);
        """
        print(text(create_index_query))
        connection.execute(text(create_index_query))


def load_seed_complements(engine, file):
    """
    example:
            engine = create_sqlalchemy_engine()
            load_seed_complements(engine, "data/tmp_cmpls_90.json.gz")
    """
    with gzip.open(file, "rt", encoding="utf-8") as f:
        #
        data = json.load(f)
    #
    df = pd.DataFrame.from_dict(data, orient="columns")
    #
    del data
    #
    # Reset index to include it as a column
    df = df.reset_index()
    #
    # Melt the DataFrame
    df_melted = df.melt(id_vars=["index"], var_name="column_name", value_name="values")
    #
    # Process values: Remove brackets and spaces
    df_melted["values"] = df_melted["values"].apply(
        lambda x: ",".join(x) if isinstance(x, list) else str(x)
    )
    #
    df_melted.columns = ["patricBeneficaryId", "patricDonorId", "seedComplements"]
    #
    df_melted_chunk = df_melted.head(2000)
    df_melted_chunk.to_sql(
        "SeedComplements", con=engine, if_exists="append", index=False
    )


def get_time():
    return time.time()


# ---------------------------


def main(database):
    # import sqlalchemy
    credentials = load_db_credentials(database)
    engine = create_sqlalchemy_engine(credentials)

    g2t = os.path.join(MAPPINGS, "update_genomeId2ncbiTaxId.tsv")
    load_tsv_to_sql(engine, g2t, "genome2taxNcbiId", if_exists="replace")

    g2p = os.path.join(MAPPINGS, "gtdb2patricIds.tsv")
    load_tsv_to_sql(engine, g2p, "patricId2genomeId", if_exists="replace")

    phen = os.path.join(ANNOTATIONS, "gtdb_phen_2023_12.tsv")
    load_tsv_to_sql(
        engine,
        phen,
        "phenDB",
        if_exists="replace",
        trim_suffix=True,
        column="gtdbId",
        escape_char="_",
    )

    # The unique_complements file does have first_column_index but we set it to False, since the table awaits for this
    upc = os.path.join(MAPPINGS, "unique_complements.tsv.gz")
    load_tsv_to_sql(
        engine, upc, "uniqueComplements", if_exists="replace", first_column_index=False
    )

    # NOTE [2025-04-22, Haris Zafeiropoulos]:
    # Pathwway complements table is ok on Luna

    # Upload pathway complementarities with complementarity ids instead of the whole thing
    compl = os.path.join(
        ANNOTATIONS, "first_10_percent.tsv.gz"
    )  # all_oneline_part_1_4.tsv.gz
    load_tsv_to_sql(
        engine,
        compl,
        "pathwayComplementarity",
        if_exists="replace",
        first_column_index=False,
    )


if __name__ == "__main__":
    # Create the parser
    parser = argparse.ArgumentParser(description="Populating the microbetag database.")

    # Add arguments
    parser.add_argument(
        "name",
        help="Name of the working environment to get credentials from .env_dev.json",
    )

    # Parse the arguments
    args = parser.parse_args()

    main(args.name)


# # import sqlalchemy
# with open(".env_dev.json", "r") as f:
#     credentials = json.load(f)
# USER_NAME, PASSWORD, HOST, DB_NAME, PORT = credentials["systems76"].values()


# engine = create_sqlalchemy_engine()

# # TODO: REPLACE WITH update_genomeId2ncbiTaxId.tsv
# g2t = os.path.join(MAPPINGS, "update_genomeId2ncbiTaxId.tsv")
# load_tsv_to_sql(engine, g2t, "genome2taxNcbiId", if_exists="replace")


# g2p = os.path.join(MAPPINGS, "gtdb2patricIds.tsv")
# load_tsv_to_sql(engine, g2p, "patricId2genomeId", if_exists="replace")


# ss = os.path.join(MAPPINGS, "update_seed_scores.tsv")
# load_seed_scores(engine, ss)

# ------------

# create_index_for_seed_scores(engine)

# load_seed_scores(engine, "scores/phylomint_scores_part_ac.gz")
