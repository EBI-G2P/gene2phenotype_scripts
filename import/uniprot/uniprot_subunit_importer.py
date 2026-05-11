#!/usr/bin/env python3

import argparse
import configparser
import csv
import re
import sys

import requests
from requests.adapters import HTTPAdapter, Retry


"""
This script imports subunit annotations from UniProtKB via the UniProt REST API.
It filters the subunit annotations to retain only those that contain keywords indicative of oligomeric state.

Usage:
    python uniprot_subunit_importer.py --config config.ini

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


URL = (
    "https://rest.uniprot.org/uniprotkb/search"
    "?query=reviewed:true+AND+organism_id:9606"
    "&fields=accession,gene_primary,xref_hgnc,cc_subunit&size=500"
)

RE_NEXT_LINK = re.compile(r'<(.+)>; rel="next"')
RETRIES = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
SESSION = requests.Session()
SESSION.mount("https://", HTTPAdapter(max_retries=RETRIES))
DEFAULT_OUTPUT_PATH = "uniprot_subunits.tsv"


def get_next_link(headers):
    if "Link" in headers:
        match = RE_NEXT_LINK.match(headers["Link"])
        if match:
            return match.group(1)
    return None


def get_batch(batch_url):
    while batch_url:
        response = SESSION.get(batch_url)
        response.raise_for_status()
        yield response
        batch_url = get_next_link(response.headers)


def extract_subunit_text(item):
    subunit_texts = []
    for comment in item.get("comments", []):
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
            else:
                continue

        filtered_part = " ".join(accepted_sentences).strip()
        if filtered_part:
            accepted_parts.append(filtered_part)
    return " | ".join(accepted_parts)


def extract_gene_symbols(item):
    genes = []
    for gene in item.get("genes", []):
        gene_symbol = gene.get("geneName", {}).get("value")
        if gene_symbol:
            genes.append(
                {
                    "gene_symbol": gene_symbol,
                    "hgnc_id": extract_hgnc_id(item, gene_symbol),
                }
            )
    return genes


def extract_hgnc_id(item, gene_symbol):
    for cross_reference in item.get("uniProtKBCrossReferences", []):
        if cross_reference.get("database") != "HGNC":
            continue

        for property_item in cross_reference.get("properties", []):
            if (
                property_item.get("key") == "GeneName"
                and property_item.get("value") == gene_symbol
            ):
                return cross_reference.get("id")

    return None


def fetch_subunit_data():
    """
    Fetches reviewed human UniProt entries with subunit annotations,
    and returns a list of dictionaries containing:
    gene symbol, HGNC ID, UniProt accession, original subunit information and filtered subunit information.
    It only returns entries where the filtered subunit information is not empty.
    """
    rows = []
    for batch in get_batch(URL):
        payload = batch.json()
        for item in payload.get("results", []):
            accession = item.get("primaryAccession")
            if not accession:
                continue

            subunit_text = extract_subunit_text(item)
            if not subunit_text:
                continue
            processed_subunit_text = filter_subunit_text(subunit_text)

            # Skip entries where the filtered subunit information is empty (i.e. no relevant keywords found)
            if not processed_subunit_text:
                continue

            for gene in extract_gene_symbols(item):
                rows.append(
                    {
                        "gene_symbol": gene["gene_symbol"],
                        "hgnc_id": gene["hgnc_id"],
                        "uniprot_accession": accession,
                        "subunit_information": subunit_text,
                        "subunit_information_filtered": processed_subunit_text,
                    }
                )
    return rows


def insert_data(g2p_db_host, g2p_db_port, g2p_db_name, g2p_user, g2p_password, subunit_data):
    for row in subunit_data:
        print("\n", row)


def main():
    parser = argparse.ArgumentParser(
        description="Export reviewed human UniProt genes with accession and subunit annotation"
    )
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

    # Fetch the subunit data from UniProt API
    subunit_data = fetch_subunit_data()

    # Insert the data into the G2P database
    insert_data(
        g2p_db_host, g2p_db_port, g2p_db_name, g2p_user, g2p_password, subunit_data
    )
    print(
        f"Exported {len(subunit_data)} gene rows with UniProt accession + subunit information.",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
