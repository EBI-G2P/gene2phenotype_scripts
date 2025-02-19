#!/usr/bin/env python3

import sys
import argparse
import MySQLdb
import configparser
from pandas import read_excel, read_csv


"""
    Script      : load_records.py

    Description : 

    Options     :

                --config : Config file containing the G2P db and API connection details (mandatory)
                        File format is the following:
                            [g2p_database]
                            host = <>
                            port = <>
                            database = <>
                            user = <> # read permission
                            password = <>
                            
                            [api]
                            api_url=<>
                
                --file : File with all records to be imported. Supported formats: csv, xlsx (mandatory)                
                --api_username: Username to connect to the G2P API (mandatory)
                --api_password: Password to connect to the G2P API (mandatory)
                --dryrun: Test script without running the import (default: 0)
"""

def read_file(file):
    file_data = None
    file_data_df = None

    # Read the file according to its format
    if file.endswith(".xlsx"):
        file_data_df = read_excel(file)
    elif file.endswith(".csv"):
        file_data_df = read_csv(file)
    else:
        sys.exit(f"ERROR: unsupported file format. Please check '{file}'")

    file_data = file_data_df.to_dict(orient='records')

    if not file_data:
        sys.exit("INFO: no data to import - input file is empty")

    # Validate the file format
    header = file_data[0].keys()
    mandatory_fields = ["gene symbol", "hgnc id", "disease name", "allelic requirement", "molecular mechanism", "confidence",
                        "publication PMID", "panel", "inferred variant consequence"]
    
    if not all(column in header for column in mandatory_fields):
        sys.exit(f"ERROR: missing data. Mandatory fields are: {mandatory_fields}")   

    return file_data

def prepare_data(data_to_load, g2p_attribs, g2p_genes, g2p_panels):
    records_to_load = {} # key: unique key identifying the record; value: json obj

    for row in data_to_load:
        valid_record = 1 # flag if record can be imported
        errors = [] # save the reason why record is invalid to print in the report

        gene_id = None
        hgnc_id = None
        disease_id = None
        genotype_id = None # allelic requirement
        molecular_mechanism_id = None
        confidence_id = None
        pmid = None
        panel_id = None
        variant_consequences = None

        # Validate the gene
        # TODO: also consider synonyms
        try:
            gene_id = g2p_genes[row["gene symbol"]]
        except KeyError:
            message = f"Gene: {row['gene symbol']} not found in G2P"
            errors.append(message)
            valid_record = 0

        # TODO: hgnc id

        # Validate the genotype (allelic requirement)
        try:
            genotype_id = g2p_attribs["genotype"][row["allelic requirement"]]
        except KeyError:
            message = f"Invalid allelic requirement (genotype): {row['allelic requirement']}"
            errors.append(message)
            valid_record = 0

        # Validate the molecular mechanism
        try:
            molecular_mechanism_id = g2p_attribs["mechanism"][row["molecular mechanism"]]
        except KeyError:
            message = f"Invalid molecular mechanism: {row['molecular mechanism']}"
            errors.append(message)
            valid_record = 0

        # Validate the confidence
        try:
            confidence_id = g2p_attribs["confidence_category"][row["confidence"]]
        except KeyError:
            message = f"Invalid confidence: {row['confidence']}"
            errors.append(message)
            valid_record = 0

        # Validate the panel
        try:
            panel_id = g2p_panels[row["panel"]]
        except KeyError:
            message = f"Invalid panel: {row['panel']}"
            errors.append(message)
            valid_record = 0

        key = row["gene symbol"]+"-"+row["disease name"]+"-"+row['allelic requirement']+"-"+row["molecular mechanism"]

        if key not in records_to_load and valid_record:
            print("\n->", row)
            

        if valid_record:
            print("Valid")
        else:
            print("Invalid:", errors)


def fetch_g2p_attribs(db_host, db_port, db_name, user, password):
    g2p_attribs = {}

    sql_attribs = """
                     SELECT a.id, a.value, at.code FROM attrib a
                     LEFT JOIN attrib_type at ON a.type_id = at.id
                     WHERE at.is_deleted = 0
                  """

    sql_mechanism = """
                        SELECT id, type, value
                        FROM cv_molecular_mechanism
                    """

    db = MySQLdb.connect(host=db_host, port=db_port, user=user, passwd=password, db=db_name)
    cursor = db.cursor()

    # Fetch all attribs
    cursor.execute(sql_attribs)
    data = cursor.fetchall()
    for row in data:
        if row[2] not in g2p_attribs:
            g2p_attribs[row[2]] = {}
            g2p_attribs[row[2]][row[1]] = row[0]
        else:
            g2p_attribs[row[2]][row[1]] = row[0]

    # Fetch all mechanisms
    cursor.execute(sql_mechanism)
    data_mechanism = cursor.fetchall()
    for row in data_mechanism:
        if row[1] not in g2p_attribs:
            g2p_attribs[row[1]] = {}
            g2p_attribs[row[1]][row[2]] = row[0]
        else:
            g2p_attribs[row[1]][row[2]] = row[0]

    return g2p_attribs

def fetch_g2p_genes(db_host, db_port, db_name, user, password):
    g2p_genes = {}

    sql_attribs = """
                     SELECT l.id, l.name FROM locus l
                     LEFT JOIN attrib a ON a.id = l.type_id
                     WHERE a.value = 'gene'
                  """

    db = MySQLdb.connect(host=db_host, port=db_port, user=user, passwd=password, db=db_name)
    cursor = db.cursor()

    # Fetch all genes
    cursor.execute(sql_attribs)
    data = cursor.fetchall()
    for row in data:
        g2p_genes[row[1]] = row[0]

    return g2p_genes

def fetch_g2p_panels(db_host, db_port, db_name, user, password):
    g2p_panels = {}

    sql_panels = """ SELECT id, name FROM panel """

    db = MySQLdb.connect(host=db_host, port=db_port, user=user, passwd=password, db=db_name)
    cursor = db.cursor()

    # Fetch all panels
    cursor.execute(sql_panels)
    data = cursor.fetchall()
    for row in data:
        g2p_panels[row[1]] = row[0]

    return g2p_panels


def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--config", required=True, help="Config file containing the API connection details")
    parser.add_argument("--file", required=True, help="File with all records to be imported")
    parser.add_argument("--api_username", required=True, help="Username to connect to the G2P API")
    parser.add_argument("--api_password", required=True, help="Password to connect to the G2P API")
    parser.add_argument("--dryrun", required=False, default=0, help="Option to test update")
    args = parser.parse_args()

    file = args.file
    config_file = args.config
    dryrun = args.dryrun
    api_username = args.api_username
    api_password = args.api_password

    # Load the config file
    config = configparser.ConfigParser()
    config.read(config_file)

    try:
        g2p_db_details = config['g2p_database']
    except KeyError:
        sys.exit("Config: G2P database details are missing from the config file")
    else:
        db_host = g2p_db_details['host']
        db_port = int(g2p_db_details['port'])
        db_name = g2p_db_details['database']
        user = g2p_db_details['user']
        password = g2p_db_details['password']

    try:
        api_url = config['api']['api_url']
    except KeyError:
        sys.exit("Config: api_url is missing from the config file")

    data_to_load = read_file(file)

    # Pre-fetch G2P attribs (genotype, mechanism, confidence, etc.)
    g2p_attribs = fetch_g2p_attribs(db_host, db_port, db_name, user, password)
    # Pre-fetch G2P genes
    g2p_genes = fetch_g2p_genes(db_host, db_port, db_name, user, password)
    # Pre-fetch G2P panels
    g2p_panels = fetch_g2p_panels(db_host, db_port, db_name, user, password)

    # Prepare the data to be imported
    prepare_data(data_to_load, g2p_attribs, g2p_genes, g2p_panels)

if __name__ == '__main__':
    main()