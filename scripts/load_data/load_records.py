#!/usr/bin/env python3

import sys
import argparse
import MySQLdb
import re
import requests
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

def prepare_data(api_url, data_to_load, g2p_attribs, g2p_genes, g2p_panels, g2p_ontology, g2p_disease_ids, g2p_mechanisms, report):
    records_to_load = {} # key: unique key identifying the record; value: json obj

    with open(report, "w") as wr:

        for row in data_to_load:
            valid_record = 1 # flag if record can be imported
            errors = [] # save the reason why record is invalid to print in the report

            disease_name = None
            pmid = None
            publication_info = None
            variant_consequences = None # column name 'inferred variant consequence'

            # Validate the gene
            # TODO: also consider synonyms
            try:
                g2p_genes[row["gene symbol"]]
            except KeyError:
                message = f"Gene: {row['gene symbol']} not found in G2P"
                errors.append(message)
                valid_record = 0

            # TODO: hgnc id

            # Validate disease name
            if f"{row['gene symbol']}-related" not in row["disease name"]:
                message = f"Invalid disease name {row['disease name']}"
                errors.append(message)
                valid_record = 0

            # Validate the genotype (allelic requirement)
            try:
                g2p_attribs["genotype"][row["allelic requirement"]]
            except KeyError:
                message = f"Invalid allelic requirement (genotype) '{row['allelic requirement']}'"
                errors.append(message)
                valid_record = 0

            # Validate the molecular mechanism
            molecular_mechanism = {}
            support = "inferred"
            try:
                g2p_mechanisms["mechanism"][None][row["molecular mechanism"]]
            except KeyError:
                message = f"Invalid molecular mechanism '{row['molecular mechanism']}'"
                errors.append(message)
                valid_record = 0
            else:
                if row["molecular mechanism evidence (by PMID)"]:
                    support = "evidence"
                molecular_mechanism = {
                    "name": row["molecular mechanism"],
                    "support": support
                }

            # Validate the mechanism synopsis
            mechanism_synopsis = {}
            if row["molecular mechanism categorisation"]:
                try:
                    g2p_mechanisms["mechanism_synopsis"][None][row["molecular mechanism categorisation"]]
                except KeyError:
                    message = f"Invalid molecular mechanism categorisation '{row['molecular mechanism categorisation']}'"
                    errors.append(message)
                    valid_record = 0
                else:
                    mechanism_synopsis = {
                        "name": row["molecular mechanism categorisation"],
                        "support": support # TODO: change
                    }

            # Validate the confidence
            try:
                g2p_attribs["confidence_category"][row["confidence"]]
            except KeyError:
                message = f"Invalid confidence '{row['confidence']}'"
                errors.append(message)
                valid_record = 0

            # Validate the panel
            try:
                g2p_panels[row["panel"]]
            except KeyError:
                message = f"Invalid panel '{row['panel']}'"
                errors.append(message)
                valid_record = 0

            # Validate the variant consequences
            # In G2P the consequence values do not include '_'
            variant_consequences = []
            for var_cons in re.split(r'\;\s*', row["inferred variant consequence"].replace("_", " ")):
                try:
                    g2p_ontology[var_cons]
                except KeyError:
                    message = f"Invalid variant consequence '{var_cons}'"
                    errors.append(message)
                    valid_record = 0
                else:
                    variant_consequences.append({
                        "support": "inferred",
                        "variant_consequence": var_cons
                    })

            # Validate the PMID - there is one publication per row
            publications = []
            if isinstance(row["publication PMID"], int):
                publication_info = fetch_pmid(api_url, row["publication PMID"])
                if not publication_info:
                    message = f"Invalid PMID {row['publication PMID']}"
                    errors.append(message)
                    valid_record = 0
                else:
                    publications.append({
                        "pmid": row["publication PMID"],
                        "year": publication_info["results"][0]["year"],
                        "title": publication_info["results"][0]["title"],
                        "source": publication_info["results"][0]["source"],
                        "authors": publication_info["results"][0]["authors"]
                    })
            else:
                message = f"Invalid PMID {row['publication PMID']}"
                errors.append(message)
                valid_record = 0

            if valid_record:
                # Append the publication comment only if the PMID is valid
                if "publication comment" in row and row["publication comment"] != "":
                    publications[0]["comment"] = row["publication comment"]

                # Append the publication number of families only if the PMID is valid
                if "publication families" in row and row["publication families"] != "":
                    publications[0]["families"] = row["publication families"]

                # Append the publication number of affected individuals only if the PMID is valid
                if "publication affected individuals" in row and row["publication affected individuals"] != "":
                    publications[0]["affectedIndividuals"] = row["publication affected individuals"]
                
                # Append the publication consanguineous only if the PMID is valid
                if "publication consanguineous" in row and row["publication consanguineous"] != "":
                    publications[0]["consanguineous"] = row["publication consanguineous"]

            # Validate the mechanism evidence (if applicable)
            # The evidence is linked to the publication
            mechanism_evidence = []
            if valid_record and "molecular mechanism evidence (by PMID)" in row and str(row["molecular mechanism evidence (by PMID)"]) != "nan":
                for m_evidence in re.split(r'\;\s*', row["molecular mechanism evidence (by PMID)"]):
                    evidence_type, evidence_value = re.split(r'\:\s*', m_evidence)
                    try:
                        g2p_mechanisms["evidence"][evidence_type.lower()][evidence_value.lower()]
                    except KeyError:
                        message = f"Invalid molecular mechanism evidence '{evidence_value}'"
                        errors.append(message)
                        valid_record = 0
                    else:
                        mechanism_evidence.append({
                            "pmid": row['publication PMID'],
                            "description": "",
                            "evidence_types": []
                        })

            # Get disease OMIM/Mondo IDs
            cross_references = []
            if row["disease mim"]:
                try:
                    omim_disease_name = g2p_disease_ids[str(row["disease mim"])]
                except KeyError:
                    message = f"Invalid OMIM ID {row['disease mim']}"
                    errors.append(message)
                    valid_record = 0
                else:
                    cross_references.append({
                        "source": "OMIM",
                        "identifier": row["disease mim"],
                        "disease_name": omim_disease_name,
                        "original_disease_name": omim_disease_name
                    })
            if row["disease MONDO"]:
                try:
                    mondo_disease_name = g2p_disease_ids[row["disease MONDO"]]
                except KeyError:
                    message = f"Invalid MONDO {row['disease MONDO']}"
                    errors.append(message)
                    valid_record = 0
                else:
                    cross_references.append({
                        "source": "Mondo",
                        "identifier": row["disease MONDO"],
                        "disease_name": mondo_disease_name,
                        "original_disease_name": mondo_disease_name
                    })

            # Get cross cutting modifier
            cross_cutting_modifier = []
            for ccm in re.split(r'\;\s*', row["cross cutting modifier"]):
                try:
                    g2p_attribs["cross_cutting_modifier"][ccm]
                except KeyError:
                    message = f"Invalid cross cutting modifier '{ccm}'"
                    errors.append(message)
                    valid_record = 0
                else:
                    cross_cutting_modifier.append(ccm)

            # Validate variant types
            variant_types = []

            # Validate phenotypes
            phenotypes = []

            # Validate variant descriptions (HGVS)
            variant_descriptions = []

            key = row["gene symbol"]+"-"+row["disease name"]+"-"+row["allelic requirement"]+"-"+row["molecular mechanism"]

            if key not in records_to_load and valid_record:
                public_comment = ""
                if row["public comment"] and str(row["public comment"]) != "nan":
                    public_comment = re.sub(r'[\n\r]+', ' ', row["public comment"])
                private_comment = ""
                if row["private comment"] and str(row["private comment"]) != "nan":
                    private_comment = re.sub(r'[\n\r]+', ' ', row["private comment"])

                # Prepare data to be send to the API
                new_record = {
                    "session_name": "",
                    "locus": row["gene symbol"],
                    "confidence": row["confidence"],
                    "panels": [row["panel"]],
                    "disease": {
                        "disease_name": row["disease name"],
                        "cross_references": cross_references
                    },
                    "allelic_requirement": row["allelic requirement"],
                    "molecular_mechanism": molecular_mechanism,
                    "cross_cutting_modifier": cross_cutting_modifier,
                    "variant_consequences": variant_consequences,
                    "phenotypes": phenotypes,
                    "publications": publications,
                    "variant_types": variant_types,
                    "public_comment": public_comment,
                    "private_comment": private_comment,
                    "mechanism_evidence": mechanism_evidence,
                    "mechanism_synopsis": mechanism_synopsis,
                    "variant_descriptions": variant_descriptions
                }

                print("To import:", new_record)

            if valid_record:
                print("Valid")
            else:
                for key in row:
                    wr.write(f"{key}: {row[key]}\t")
                wr.write(f"Cannot load record: {';'.join(errors)}\n")


def fetch_g2p_attribs(db_host, db_port, db_name, user, password):
    g2p_attribs = {}

    sql_attribs = """
                     SELECT a.id, a.value, at.code FROM attrib a
                     LEFT JOIN attrib_type at ON a.type_id = at.id
                     WHERE at.is_deleted = 0
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

    return g2p_attribs

def fetch_g2p_mechanisms(db_host, db_port, db_name, user, password):
    g2p_mechanisms = {}

    sql_mechanism = """
                        SELECT id, type, value, subtype
                        FROM cv_molecular_mechanism
                    """

    db = MySQLdb.connect(host=db_host, port=db_port, user=user, passwd=password, db=db_name)
    cursor = db.cursor()

    # Fetch all mechanisms
    cursor.execute(sql_mechanism)
    data_mechanism = cursor.fetchall()
    for row in data_mechanism:
        if row[1] not in g2p_mechanisms:
            g2p_mechanisms[row[1]] = {}
            g2p_mechanisms[row[1]][row[3]] = {}
            g2p_mechanisms[row[1]][row[3]][row[2]] = row[0]
        elif row[3] not in g2p_mechanisms[row[1]]:
            g2p_mechanisms[row[1]][row[3]] = {}
            g2p_mechanisms[row[1]][row[3]][row[2]] = row[0]
        else:
            g2p_mechanisms[row[1]][row[3]][row[2]] = row[0]

    return g2p_mechanisms

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

def fetch_g2p_ontology(db_host, db_port, db_name, user, password):
    g2p_ontology = {}

    sql_ontology = """ SELECT o.id, o.term FROM ontology_term o
                       LEFT JOIN attrib a ON a.id = o.group_type_id
                       WHERE a.value = 'variant_type'
                   """

    db = MySQLdb.connect(host=db_host, port=db_port, user=user, passwd=password, db=db_name)
    cursor = db.cursor()

    # Fetch all panels
    cursor.execute(sql_ontology)
    data = cursor.fetchall()
    for row in data:
        g2p_ontology[row[1]] = row[0]

    return g2p_ontology

def fetch_g2p_disease_ids(db_host, db_port, db_name, user, password):
    g2p_ext_ids = {}

    sql_disease = """ SELECT identifier, disease
                      FROM gene_disease
                  """

    db = MySQLdb.connect(host=db_host, port=db_port, user=user, passwd=password, db=db_name)
    cursor = db.cursor()

    # Fetch all panels
    cursor.execute(sql_disease)
    data = cursor.fetchall()
    for row in data:
        g2p_ext_ids[row[0]] = row[1]

    return g2p_ext_ids

def fetch_pmid(api_url, pmid):
    result = None
    publication_url = api_url+"/publication/"+str(pmid)

    response = requests.get(publication_url)
    if response.status_code == 200:
        result = response.json()

    return result


def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--config", required=True, help="Config file containing the API connection details")
    parser.add_argument("--file", required=True, help="File with all records to be imported")
    parser.add_argument("--api_username", required=True, help="Username to connect to the G2P API")
    parser.add_argument("--api_password", required=True, help="Password to connect to the G2P API")
    parser.add_argument("--report", required=False, default="report.txt", help="Report file")
    parser.add_argument("--dryrun", required=False, default=0, help="Option to test update")
    args = parser.parse_args()

    file = args.file
    config_file = args.config
    report = args.report
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
    # Pre-fetch G2P ontology terms
    g2p_ontology = fetch_g2p_ontology(db_host, db_port, db_name, user, password)
    # Pre-fetch G2P OMIM/Mondo IDs
    g2p_disease_ids = fetch_g2p_disease_ids(db_host, db_port, db_name, user, password)
    # Pre-fetch G2P mechanism attribs
    g2p_mechanisms = fetch_g2p_mechanisms(db_host, db_port, db_name, user, password)

    # Prepare the data to be imported
    prepare_data(api_url, data_to_load, g2p_attribs, g2p_genes, g2p_panels, g2p_ontology, g2p_disease_ids, g2p_mechanisms, report)

if __name__ == '__main__':
    main()