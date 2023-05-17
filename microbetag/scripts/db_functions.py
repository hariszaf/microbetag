import os, sys
import decimal
import pandas as pd
import mysql.connector
from .variables import *


def check_connection_to_db(db_user, db_password, db_host, db_database):

    try:
        db = mysql.connector.connect(user=db_user, password=db_password, host=db_host, database=db_database)
        cursor = db.cursor()        
        cursor.execute("SELECT VERSION()")
        results = cursor.fetchone()

        if results:
            return True
        else:
            return False               

    except mysql.connector.Error:
        return False


def query_to_microbetagDB(phrase):
    """
    Function to get functional traits as assigned in a genome from the phenotrex software
    phenotrex.readthedocs.io/ and stored in the microbetagDB.
    """
    try:        
        cnx = mysql.connector.connect(user=USER_NAME, password=PASSWORD, host=HOST, database=DB_NAME)
        cnx.get_warnings = True
        cursor = cnx.cursor()
        cursor.execute(phrase)
        res = []
        for row in cursor:
            res.append(row)

        warnings = cursor.fetchwarnings()
        if warnings: 
            for i in range(len(warnings)):
                print("\t ** Warning - "+warnings[i][2])
        cursor.close()
        cnx.commit()
        cnx.close()
        return res

    except mysql.connector.Error as err:
        print("Something went wrong: {}".format(err))
        print(phrase)

def get_phenDB_traits(genome_id):
    phrase = "SELECT * FROM phenDB WHERE gtdbId = '" + genome_id + "';"
    phendb_traits = list(query_to_microbetagDB(phrase)[0])
    phendb_traits = [str(x) if isinstance(x, decimal.Decimal) else x for x in phendb_traits]
    return phendb_traits

def get_column_names():
    phrase = "SELECT `COLUMN_NAME` FROM `INFORMATION_SCHEMA`.`COLUMNS` WHERE `TABLE_SCHEMA`='" +\
             DB_NAME + "' AND `TABLE_NAME`='phenDB';"
    colnames = query_to_microbetagDB(phrase)
    colnames = [x[0] for x in colnames]
    return colnames

