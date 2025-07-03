#!/usr/bin/env python3

import sys
import argparse
import configparser
import csv
import MySQLdb
import os.path
import requests
from lxml import etree as ET


"""
    Description: Script to generate XML file for EBI Search

    Params:
            --config : Config file name containing the database and API connection info (mandatory)
                    File format is the following: 
                        [g2p_database]
                        host = <>
                        port = <>
                        user = <>
                        password = <>
                        name = <>

                        [api]
                        api_url = <>

            --output_dir : Path to the output directory where XML file is going to be saved (mandatory)
    """


def get_g2p_version(api_url: str) -> str:
    """
    Fetch the latest G2P version from the API.

    Args:
        api_url (str): G2P API URL

    Returns:
        str: G2P version
    """
    g2p_version = None
    url = api_url + "/reference_data/"

    try:
        response = requests.get(url)
    except Exception as e:
        sys.exit(f"Error while fetching the G2P version:", e)
    else:
        if response.status_code == 200:
            result = response.json()
            for obj in result:
                if obj["key"] == "g2p_release":
                    g2p_version = obj["version"]
        else:
            sys.exit(f"Failed to fetch G2P version from API")

    return g2p_version


def dump_g2p_records(
    api_url: str, db_host: str, db_port: int, db_name: str, user: str, password: str
) -> dict[str, dict]:
    """
    Queries the G2P API to download all records associated with public panels.
    It also queries the database to fetch the gene Ensembl IDs as this ID is not
    available in the download file fetched from the API.

    Args:
        api_url (str): G2P API URL
        db_host (str): G2P database host name
        db_port (int): G2P database port number
        db_name (str): G2P database name
        user (str): user with read access
        password (str): password

    Returns:
        dict[str, dict]: All records associated with public panels
    """
    url = api_url + "/panel/all/download/"
    g2p_file = "g2p_data.txt"
    records = {}
    genes = {}

    # Query the database to fetch the Ensembl gene IDs
    try:
        db = MySQLdb.connect(
            host=db_host, port=db_port, user=user, passwd=password, db=db_name
        )
    except MySQLdb.Error as error:
        sys.exit("Database connection failed:", error)
    else:
        sql = """ SELECT l.name, li.identifier
                  FROM locus l
                  LEFT JOIN locus_identifier li ON li.locus_id = l.id
                  LEFT JOIN source s ON s.id = li.source_id
                  WHERE s.name = "Ensembl"
              """

        cursor = db.cursor()
        cursor.execute(sql)
        data = cursor.fetchall()
        for row in data:
            genes[row[0]] = row[1]
        cursor.close()
        db.close()

    # Query the API to download file
    try:
        response = requests.get(url, stream=True)
    except Exception as e:
        sys.exit(f"Error while downloading G2P data:", e)
    else:
        if response.status_code == 200:
            with open(g2p_file, "wb") as wr:
                for chunk in response.iter_content(chunk_size=128):
                    wr.write(chunk)
        else:
            sys.exit(f"Failed to download G2P data")

    # Read the file
    with open(g2p_file, "r") as fh:
        reader = csv.DictReader(fh)
        for line in reader:
            g2p_id = line["g2p id"]

            if g2p_id not in records:
                records[g2p_id] = {
                    "disease": line["disease name"],
                    "confidence": line["confidence"],
                    "genotype": line["allelic requirement"],
                    "mechanism": line["molecular mechanism"],
                    "gene": line["gene symbol"],
                    "panels": line["panel"],
                    "gene_ensembl": genes[line["gene symbol"]],
                    "hgnc_id": line["hgnc id"],
                    "disease_mondo": line["disease MONDO"],
                    "disease_mim": line["disease mim"],
                    "pmids": line["publications"],
                }
            else:
                print(f"WARNING: duplicated record '{g2p_id}'")

    return records


def add_element(parent: object, name: str, value: str) -> None:
    """
    Add a field to the XML obj.

    Args:
        parent (object): ET element obj
        name (str): name of the element field
        value (str): value of the field
    """
    if value:
        child = ET.SubElement(parent, name)
        child.text = value


def create_xml(g2p_version: str, g2p_records: dict[str, dict]) -> bytes:
    """
    Generate the XML data with all the records that we want to display on the
    EBI search engine.

    Args:
        g2p_version (str): G2P version
        g2p_records (dict[str, dict]): List of all records linked to public panels

    Returns:
        bytes: xml data generated with the G2P records info
    """
    database_element = ET.Element("database")

    # Add metadata
    add_element(database_element, "name", "G2P")
    add_element(
        database_element,
        "description",
        "Gene2Phenotype (G2P) is a detailed collection of expert curated gene-disease associations with information on allelic requirement, observed variant classes and disease mechanism",
    )
    add_element(database_element, "url", "https://www.ebi.ac.uk/gene2phenotype/")
    add_element(database_element, "url_search", "https://www.ebi.ac.uk/gene2phenotype/lgd/")
    add_element(database_element, "release", g2p_version)
    add_element(database_element, "entry_count", str(len(g2p_records)))

    # Add records
    entries_element = ET.SubElement(database_element, "entries")
    for entry in g2p_records:
        entry_element = ET.SubElement(entries_element, "entry", id=entry, acc=entry)
        entry_name = (
            f"{g2p_records[entry]['disease']} ({g2p_records[entry]['confidence']})"
        )
        add_element(entry_element, "name", entry_name)
        # Additional fields
        add_fields_element = ET.SubElement(entry_element, "additional_fields")
        f = ET.SubElement(add_fields_element, "field", name="gene")
        f.text = g2p_records[entry]["gene"]
        f = ET.SubElement(add_fields_element, "field", name="disease")
        f.text = g2p_records[entry]["disease"]
        f = ET.SubElement(add_fields_element, "field", name="genotype")
        f.text = g2p_records[entry]["genotype"]
        f = ET.SubElement(add_fields_element, "field", name="mechanism")
        f.text = g2p_records[entry]["mechanism"]
        f = ET.SubElement(add_fields_element, "field", name="confidence")
        f.text = g2p_records[entry]["confidence"]

        # Cross references - gene ID
        xrefs_element = ET.SubElement(entry_element, "cross_references")
        # Gene IDs: Ensembl and HGNC
        ET.SubElement(
            xrefs_element, "ref", dbname="HGNC", dbkey=g2p_records[entry]["hgnc_id"]
        )
        ET.SubElement(
            xrefs_element,
            "ref",
            dbname="ENSEMBL_GENE",
            dbkey=g2p_records[entry]["gene_ensembl"],
        )

        # Cross references - disease ontology
        if g2p_records[entry]["disease_mondo"] != "":
            ET.SubElement(
                xrefs_element,
                "ref",
                dbname="Mondo",
                dbkey=g2p_records[entry]["disease_mondo"],
            )
        if g2p_records[entry]["disease_mim"] != "":
            ET.SubElement(
                xrefs_element,
                "ref",
                dbname="OMIM_DISEASE",
                dbkey=g2p_records[entry]["disease_mim"],
            )

        # Cross references - publications
        if g2p_records[entry]["pmids"] != "":
            publications_list = g2p_records[entry]["pmids"].split("; ")
            for pmid in publications_list:
                ET.SubElement(xrefs_element, "ref", dbname="EUROPEPMC", dbkey=pmid)

    return ET.tostring(
        database_element, pretty_print=True, encoding="UTF-8", xml_declaration=True
    )


def main():
    parser = argparse.ArgumentParser(
        description="Generate a G2P XML file for the EBI search engine"
    )
    parser.add_argument(
        "--config",
        required=True,
        help="Config file containing details to the G2P database",
    )
    parser.add_argument(
        "--output_dir", required=True, help="Path to the output directory"
    )
    args = parser.parse_args()

    config_file = args.config
    output_dir = args.output_dir

    if not os.path.exists(output_dir):
        sys.exit(f"Invalid output directory '{output_dir}'")

    # Load the config file
    config = configparser.ConfigParser()
    config.read(config_file)

    try:
        g2p_config = config["g2p_database"]
    except KeyError:
        sys.exit("ERROR: 'g2p_database' missing from config file")
    else:
        db_host = g2p_config["host"]
        db_port = g2p_config["port"]
        db_name = g2p_config["name"]
        user = g2p_config["user"]
        password = g2p_config["password"]

    try:
        api = config["api"]
    except KeyError:
        sys.exit("ERROR: 'api' missing from config file")
    else:
        api_url = api["api_url"]

    # Get the G2P version from the meta table
    g2p_version = get_g2p_version(api_url)
    # Fetch all records to be available in the EBI search
    g2p_records = dump_g2p_records(
        api_url, db_host, int(db_port), db_name, user, password
    )
    # # Generate the XML with all records
    xml_data = create_xml(g2p_version, g2p_records)

    # Write to the output file - open in binary mode because 'xml_data' is a bytes obj
    output_file = os.path.join(output_dir, "g2p_records.xml")

    with open(output_file, "wb") as wr:
        wr.write(xml_data)


if __name__ == "__main__":
    main()
