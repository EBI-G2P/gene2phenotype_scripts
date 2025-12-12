#!/usr/bin/env python3

import argparse
import configparser
import datetime
import re
import sys

import MySQLdb
import requests
from requests.adapters import HTTPAdapter, Retry


"""
This script imports required gene information from Uniprot via the Uniprot REST API.

Usage: python uniprot_importer.py [OPTIONS]

Options:
        --config    Config file with details to the G2P database (mandatory)
                        File format is the following:
                                [g2p_database]
                                host = <>
                                port = <>
                                user = <>
                                password = <>
                                name = <>
"""


# Uniprot data fetch URL
url = "https://rest.uniprot.org/uniprotkb/search?query=reviewed:true+AND+organism_id:9606&fields=accession,cc_function,xref_mim,xref_hgnc,gene_primary&size=500"

# Configuration to fetch Uniprot data
re_next_link = re.compile(r'<(.+)>; rel="next"')
retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))

# Global variable to store Uniprot release version
uniprot_release = None


def get_next_link(headers):
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)


def get_batch(batch_url):
    global uniprot_release
    while batch_url:
        response = session.get(batch_url)
        response.raise_for_status()
        if uniprot_release is None:
            uniprot_release = response.headers["X-UniProt-Release"]
        yield response
        batch_url = get_next_link(response.headers)


def get_database_value(database, dataItem, gene_symbol):
    if "uniProtKBCrossReferences" in dataItem and len(
        dataItem["uniProtKBCrossReferences"]
    ):
        for item in dataItem["uniProtKBCrossReferences"]:
            # Example: {'database': 'HGNC', 'id': 'HGNC:4764', 'properties': [{'key': 'GeneName', 'value': 'H3-3A'}]
            if database == "HGNC" and item["properties"][0]["value"] == gene_symbol:
                return item["id"]

            # Example: {'database': 'MIM', 'id': '601058', 'properties': [{'key': 'Type', 'value': 'gene'}]}
            elif database == "MIM" and item["properties"][0]["value"] == "gene":
                return item["id"]

            # This checks if there is a HGNC ID
            elif item["database"] == database and gene_symbol is None:
                return item["id"]

    return None


def is_protein_function_available(dataItem):
    return (
        "comments" in dataItem
        and len(dataItem["comments"]) != 0
        and "texts" in dataItem["comments"][0]
        and len(dataItem["comments"][0]["texts"]) != 0
        and "value" in dataItem["comments"][0]["texts"][0]
    )


# Function to fetch Uniprot data
def fetch_all_data():
    total_items = []
    for batch in get_batch(url):
        current_batch_json = batch.json()
        for item in current_batch_json["results"]:
            # If protein function or HGNC id is not available then don't consider the data entry
            if (
                is_protein_function_available(item)
                and get_database_value("HGNC", item, None) is not None
            ):
                # 'genes' is a list which can have multiple values (gene symbols)
                for gene in item["genes"]:
                    current_item = {}
                    current_item["gene_symbol"] = gene["geneName"]["value"]
                    current_item["accession"] = item["primaryAccession"]
                    current_item["protein_function"] = item["comments"][0]["texts"][0][
                        "value"
                    ]
                    current_item["HGNC"] = get_database_value(
                        "HGNC", item, current_item["gene_symbol"]
                    )
                    current_item["MIM"] = get_database_value(
                        "MIM", item, current_item["gene_symbol"]
                    )
                    total_items.append(current_item)

    print(
        f"Uniprot data successfully fetched via Uniprot API ({len(total_items)} entries)"
    )
    return total_items


# Function to insert Uniprot data to database
def insert_uniprot_data(db_host, db_port, db_name, db_user, db_password, total_items):
    sql_truncate = """ TRUNCATE TABLE uniprot_annotation """

    sql_count = """ SELECT COUNT(*) from uniprot_annotation """

    sql_source = (
        """ SELECT id, name FROM source WHERE name = 'UniProt' or name = 'HGNC' """
    )

    sql_identifier = (
        """ SELECT identifier, locus_id FROM locus_identifier WHERE source_id = %s"""
    )

    insert_sql = """ INSERT INTO uniprot_annotation(uniprot_accession, gene_id, hgnc, gene_symbol, mim, protein_function, source_id)
                  VALUES(%s, %s, %s, %s, %s, %s, %s) """

    sql_meta = """ INSERT INTO meta(`key`, date_update, is_public, description, source_id, version)
                    VALUES(%s,%s,%s,%s,%s,%s) """

    db = MySQLdb.connect(
        host=db_host, port=db_port, user=db_user, passwd=db_password, db=db_name
    )
    cursor = db.cursor()

    # Save number of rows from uniprot_annotation
    cursor.execute(sql_count)
    previous_number_rows = cursor.fetchone()

    # Truncate uniprot_annotation table before the import is run
    cursor.execute(sql_truncate)
    db.commit()

    # Fetch source ids
    source_ids = {}
    cursor.execute(sql_source)
    data = cursor.fetchall()
    if len(data) != 0:
        for row in data:
            source_ids[row[1]] = row[0]
    # Fetch locus identifiers
    identifier_to_locus_id_map = {}
    cursor.execute(sql_identifier, [source_ids["HGNC"]])
    data = cursor.fetchall()
    if len(data) != 0:
        for row in data:
            identifier_to_locus_id_map[row[0]] = row[1]
    # Insert Uniprot data
    insert_count = 0
    for item in total_items:
        if item["HGNC"] in identifier_to_locus_id_map:
            cursor.execute(
                insert_sql,
                [
                    item["accession"],
                    identifier_to_locus_id_map[item["HGNC"]],
                    item["HGNC"],
                    item["gene_symbol"],
                    item["MIM"],
                    item["protein_function"],
                    source_ids["UniProt"],
                ],
            )
            insert_count += 1
    # Insert import info into meta table
    cursor.execute(
        sql_meta,
        [
            "import_uniprot",
            datetime.datetime.now(),
            0,
            "Import Uniprot data",
            source_ids["UniProt"],
            uniprot_release,
        ],
    )
    db.commit()
    db.close()
    print("Uniprot data successfully inserted into G2P database.")
    print(f"Previous total number of uniprot entries: {previous_number_rows[0]}")
    print(f"Total Uniprot entries fetched: {len(total_items)}")
    print(f"Total Uniprot entries inserted: {insert_count}")
    print(
        "Note: Only Uniprot data entries with existing Gene information in the database will be inserted."
    )


def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument(
        "--config", required=True, help="Config file with details to the G2P database"
    )

    args = parser.parse_args()
    config_file = args.config

    # Load the config file
    config = configparser.ConfigParser()
    config.read(config_file)

    try:
        g2p_config = config["g2p_database"]
    except KeyError:
        sys.exit("ERROR: 'g2p_database' missing from config file")
    else:
        g2p_db_host = g2p_config["host"]
        g2p_db_port = int(g2p_config["port"])
        g2p_db_name = g2p_config["name"]
        g2p_user = g2p_config["user"]
        g2p_password = g2p_config["password"]

    total_items = fetch_all_data()
    insert_uniprot_data(
        g2p_db_host, g2p_db_port, g2p_db_name, g2p_user, g2p_password, total_items
    )


if __name__ == "__main__":
    main()
