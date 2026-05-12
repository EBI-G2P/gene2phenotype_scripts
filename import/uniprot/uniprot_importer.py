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
Data retrieved from Uniprot includes:
 - protein function
 - subunit information with filtering to keep only relevant information about the quaternary structure (e.g. homodimer, heterodimer, monomer, etc.)

Usage: python uniprot_importer.py --config <config_file>

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
url = "https://rest.uniprot.org/uniprotkb/search?query=reviewed:true+AND+organism_id:9606&fields=accession,cc_function,cc_subunit,xref_hgnc,gene_primary&size=500"

# Configuration to fetch Uniprot data
re_next_link = re.compile(r'<(.+)>; rel="next"')
retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))

# Global variable to store Uniprot release version
uniprot_release = None


def fetch_uniprot_release():
    global uniprot_release
    release_url = "https://rest.uniprot.org/uniprotkb/search?query=reviewed:true+AND+organism_id:9606&fields=accession&size=1"
    response = session.get(release_url)
    response.raise_for_status()
    uniprot_release = response.headers["X-UniProt-Release"]
    return uniprot_release


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


def get_hgnc_id(dataItem, gene_symbol=None):
    if "uniProtKBCrossReferences" in dataItem and len(
        dataItem["uniProtKBCrossReferences"]
    ):
        for item in dataItem["uniProtKBCrossReferences"]:
            if item["database"] != "HGNC":
                continue

            # This checks if there is a HGNC ID
            if gene_symbol is None:
                return item["id"]

            # Example: {'database': 'HGNC', 'id': 'HGNC:4764', 'properties': [{'key': 'GeneName', 'value': 'H3-3A'}]
            if any(
                property_item["key"] == "GeneName"
                and property_item["value"] == gene_symbol
                for property_item in item.get("properties", [])
            ):
                return item["id"]

    return None


def extract_protein_function_text(dataItem):
    for comment in dataItem.get("comments", []):
        if comment.get("commentType") != "FUNCTION":
            continue
        for text_item in comment.get("texts", []):
            value = text_item.get("value")
            if value:
                return value.replace("\n", " ").strip()
    return None


def extract_subunit_text(dataItem):
    subunit_texts = []
    for comment in dataItem.get("comments", []):
        if comment.get("commentType") != "SUBUNIT":
            continue
        for text_item in comment.get("texts", []):
            value = text_item.get("value")
            if value:
                subunit_texts.append(value.replace("\n", " ").strip())
    return " | ".join(subunit_texts)


def filter_subunit_text(subunit_text):
    keywords = [
        "homodimer",
        "heterodimer",
        "monomer",
        "dimer",
        "pentamer",
        "pentameric",
        "hexamer",
        "hexameric",
        "octamer",
        "dodecamer",
        "hexadecamer",
        "homomer",
        "homodimer",
        "homotrimer",
        "homotetramer",
        "homopentamer",
        "homohexamer",
        "homoheptamer",
        "homodecamer",
        "homododecamer",
        "homomultimer",
        "heterodimer",
        "heterotrimer",
        "heterotetramer",
        "heteropentamer",
        "heterohexamer",
        "heteromultimer",
        "heteromultimeric",
        "homooligomer",
        "homooligomerizes",
        "oligomer",
        "oligomerize",
        "self-associates",
    ]

    if not subunit_text:
        return ""

    accepted_parts = []
    for part in subunit_text.split(" | "):
        part = part.strip()
        if not part:
            continue

        accepted_sentences = []
        for sentence in re.split(r"(?<=\.)\s+", part):
            sentence = sentence.strip()
            if not sentence:
                continue
            if any(s in sentence for s in keywords):
                accepted_sentences.append(sentence)

        filtered_part = " ".join(accepted_sentences).strip()
        if filtered_part:
            accepted_parts.append(filtered_part)
    return " | ".join(accepted_parts)


# Function to fetch Uniprot data
def fetch_all_data():
    total_items = []
    for batch in get_batch(url):
        current_batch_json = batch.json()
        for item in current_batch_json["results"]:
            # If HGNC id is not available then don't consider the data entry.
            if get_hgnc_id(item) is None:
                continue

            accession = item["primaryAccession"]
            protein_function_text = extract_protein_function_text(item)
            filtered_subunit_text = filter_subunit_text(extract_subunit_text(item))

            if not protein_function_text and not filtered_subunit_text:
                continue

            # 'genes' is a list which can have multiple values (gene symbols).
            for gene in item["genes"]:
                gene_symbol = gene["geneName"]["value"]
                hgnc_id = get_hgnc_id(item, gene_symbol)

                if protein_function_text:
                    total_items.append(
                        {
                            "accession": accession,
                            "HGNC": hgnc_id,
                            "annotation": protein_function_text,
                            "data_type": "uniprot_protein_function",
                        }
                    )

                if filtered_subunit_text:
                    total_items.append(
                        {
                            "accession": accession,
                            "HGNC": hgnc_id,
                            "annotation": filtered_subunit_text,
                            "data_type": "uniprot_subunit",
                        }
                    )

    print(
        f"Uniprot data successfully fetched via Uniprot API ({len(total_items)} entries)"
    )
    return total_items


def get_current_uniprot_release(db_host, db_port, db_name, db_user, db_password):
    sql_meta = """ SELECT m.version
                   FROM meta m
                   JOIN source s ON s.id = m.source_id
                   WHERE m.`key` = %s AND s.name = %s
                   ORDER BY m.date_update DESC
                   LIMIT 1 """

    db = MySQLdb.connect(
        host=db_host, port=db_port, user=db_user, passwd=db_password, db=db_name
    )
    cursor = db.cursor()
    cursor.execute(sql_meta, ["import_uniprot", "UniProt"])
    data = cursor.fetchone()
    db.close()

    if data is None:
        return None

    return data[0]


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

    sql_attrib = """ SELECT a.id
                     FROM attrib a
                     JOIN attrib_type at ON at.id = a.type_id
                     WHERE a.value = %s AND at.code = %s """

    insert_sql = """ INSERT INTO uniprot_annotation(uniprot_accession, gene_id, annotation, source_id, data_type_id)
                  VALUES(%s, %s, %s, %s, %s) """

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

    data_type_ids = {}
    for data_type in ["uniprot_protein_function", "uniprot_subunit"]:
        cursor.execute(sql_attrib, [data_type, "uniprot_data"])
        data_type_row = cursor.fetchone()
        if data_type_row is None:
            sys.exit(f"ERROR: attrib '{data_type}' of type 'uniprot_data' is missing")
        data_type_ids[data_type] = data_type_row[0]

    # Fetch locus identifiers
    identifier_to_locus_id_map = {}
    cursor.execute(sql_identifier, [source_ids["HGNC"]])
    data = cursor.fetchall()
    if len(data) != 0:
        for row in data:
            identifier_to_locus_id_map[row[0]] = row[1]
    # Insert Uniprot data
    insert_count = 0
    insert_count_by_attrib = {}
    for item in total_items:
        if item["HGNC"] in identifier_to_locus_id_map:
            cursor.execute(
                insert_sql,
                [
                    item["accession"],
                    identifier_to_locus_id_map[item["HGNC"]],
                    item["annotation"],
                    source_ids["UniProt"],
                    data_type_ids[item["data_type"]],
                ],
            )
            insert_count += 1
            insert_count_by_attrib[item["data_type"]] = (
                insert_count_by_attrib.get(item["data_type"], 0) + 1
            )
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
    print("Total Uniprot entries inserted by attrib:")
    for data_type in sorted(insert_count_by_attrib):
        print(f"  {data_type}: {insert_count_by_attrib[data_type]}")
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

    latest_uniprot_release = fetch_uniprot_release()
    current_uniprot_release = get_current_uniprot_release(
        g2p_db_host, g2p_db_port, g2p_db_name, g2p_user, g2p_password
    )

    if current_uniprot_release == latest_uniprot_release:
        print(
            f"Skipping update: UniProt data is up-to-date ({latest_uniprot_release})."
        )
        return

    total_items = fetch_all_data()
    insert_uniprot_data(
        g2p_db_host, g2p_db_port, g2p_db_name, g2p_user, g2p_password, total_items
    )


if __name__ == "__main__":
    main()
