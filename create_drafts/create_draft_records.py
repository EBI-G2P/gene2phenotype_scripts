#!/usr/bin/env python3

import sys
import argparse
import configparser
import csv
import io
import re
import os.path
import requests
import json
from typing import List, Dict


"""
Script to read the pre-draft data from ClinGen or PanelApp JSON file and create draft records.

Options:
        --source         ClinGen or PanelApp
        --input_file     Data to be used to create drafts (supported format: json)

Example how to run:
    python create_draft_records.py --source ClinGen --input_file clingen_reviewed_data.json
"""

valid_mechanisms = [
    "loss of function",
    "gain of function",
    "dominant negative",
    "undetermined non loss of function",
    "undetermined",
]


def login(
    api_username: str, api_password: str, api_url: str
) -> requests.cookies.RequestsCookieJar:
    """Login into G2P API"""
    login_url = f"{api_url.rstrip('/')}/login/"

    response = requests.post(
        login_url, json={"username": api_username, "password": api_password}
    )

    if response.status_code != 200:
        sys.exit("Login failed. Check your credentials and API URL.")

    return response.cookies


def logout(api_url: str, cookies: requests.cookies.RequestsCookieJar) -> None:
    """Logout of the API"""
    logout_url = f"{api_url.rstrip('/')}/logout/"

    response = requests.post(logout_url, cookies=cookies)

    if response.status_code != 204:
        sys.exit("Logout failed. Check your credentials and API URL.")


def get_all_automatic_drafts(
    api_url: str, cookies: requests.cookies.RequestsCookieJar
) -> dict:
    """Fetch all automatic draft records from G2P"""
    drafts = {}

    url = f"{api_url.rstrip('/')}/curations/?scope=all&type=automatic"
    response = requests.get(url, cookies=cookies)
    if response.status_code != 200:
        sys.exit(f"Failed to fetch draft records: {response.text}")
    else:
        data_response = response.json()["results"]
        for draft in data_response:
            drafts[draft["session_name"]] = draft

    return drafts


def insert_draft_record(
    api_url: str, cookies: requests.cookies.RequestsCookieJar, draft_data: dict
) -> None:
    """
    Insert a draft record into G2P
    The endpoint accepts the following data format:
        {
            "json_data": {...},
            "status": "automatic"  # optional, default is 'manual'
        }
    """
    url = f"{api_url.rstrip('/')}/add/curation/"
    response = requests.post(url, json=draft_data, cookies=cookies)
    if response.status_code != 200:
        sys.exit(f"Failed to insert draft record: {response.text}")
    else:
        print(f"Draft record inserted successfully: {response.json()}")


def fetch_pmid_info(api_url: str, pmid: str) -> dict:
    url = f"{api_url.rstrip('/')}/publication/{pmid}"
    response = requests.get(url)
    if response.status_code != 200:
        sys.exit(f"Cannot fetch PMID {pmid}")
    else:
        return response.json()["results"][0]


def fetch_disease_info(api_url: str, disease_id: str) -> dict:
    url = f"{api_url.rstrip('/')}/external_disease/{disease_id}"
    response = requests.get(url)
    if response.status_code != 200:
        sys.exit(f"Cannot fetch disease {disease_id}")
    else:
        return response.json()["results"][0]


def prepare_draft_records(
    input_file: str,
    source: str,
    automatic_drafts: dict,
    api_url: str,
    cookies: requests.cookies.RequestsCookieJar,
) -> None:
    """
    Read the ClinGen json file
    """
    existing_records = []

    with open(input_file) as fh:
        data = json.load(fh)
        for record in data:
            final_draft = {}

            if record["status"] == "approved" or record["status"] == "needs changes":
                if "g2p_data" in record["comments"] and (
                    "matches the existing g2p" in record["comments"]["g2p_data"].lower()
                    or "draft matches" in record["comments"]["g2p_data"].lower()
                ):
                    existing_records.append(record)
                else:
                    final_draft["extra_comment"] = ""
                    # Add an extra comment related to the necessary changes
                    if record["status"] == "needs changes":
                        final_draft["extra_comment"] += (
                            "Note: this draft needs changes" + "\n"
                        )
                    if "g2p_data" in record["comments"]:
                        final_draft["extra_comment"] += (
                            "Existing G2P data: "
                            + record["comments"]["g2p_data"]
                            + "\n"
                        )

                    final_draft["locus"] = record["gene_symbol"]
                    final_draft["panels"] = ["Ear disorders"]
                    final_draft["public_comment"] = ""
                    final_draft["private_comment"] = ""
                    final_draft["cross_cutting_modifier"] = []
                    final_draft["variant_types"] = []
                    final_draft["variant_descriptions"] = []
                    final_draft["variant_consequences"] = []

                    # Add source details to JSON
                    final_draft["source_data"] = {"name": source}
                    if "url" in record:
                        final_draft["source_data"]["url"] = record["url"]

                    final_draft["confidence"] = ""
                    # Save the confidence in the source field
                    final_draft["source_data"]["confidence"] = record["confidence"].lower()
                    if "confidence" in record["comments"]:
                        final_draft["extra_comment"] += (
                            "Confidence comment: "
                            + record["comments"]["confidence"]
                            + "\n"
                        )

                    final_draft["publications"] = []
                    for pmid in record["pmids"]:
                        publication_data = fetch_pmid_info(api_url, pmid)
                        final_draft["publications"].append(
                            {
                                "pmid": pmid,
                                "year": publication_data["year"],
                                "title": publication_data["title"],
                                "source": publication_data["source"],
                                "authors": publication_data["authors"],
                                "comment": "",
                                "families": None,
                                "ancestries": "",
                                "consanguineous": "unknown",
                                "affectedIndividuals": None,
                            }
                        )

                    final_draft["phenotypes"] = []
                    if record["phenotypes"]:
                        final_draft["extra_comment"] += (
                            "Phenotypes: " + "; ".join(record["phenotypes"]) + "\n"
                        )

                    # Format allelic requirement
                    final_draft["allelic_requirement"] = ""
                    if record["allelic_requirement"] == "autosomal recessive":
                        final_draft["allelic_requirement"] = "biallelic_autosomal"
                    elif record["allelic_requirement"] == "autosomal dominant":
                        final_draft["allelic_requirement"] = "monoallelic_autosomal"
                    elif record["allelic_requirement"] == "X-linked":
                        final_draft["allelic_requirement"] = "monoallelic_X"
                    else:
                        final_draft["extra_comment"] += (
                            "Unsupported allelic requirement: "
                            + record["allelic_requirement"]
                            + "\n"
                        )

                    # Mechanism
                    final_draft["source_data"]["mechanism_comment"] = ""
                    # Check if the mechanism value is valid, if not add it to the extra_comment
                    valid_mechanism = ""
                    if (
                        record["mechanism"] != ""
                        and record["mechanism"].lower().replace("-", " ")
                        not in valid_mechanisms
                    ):
                        final_draft["source_data"]["mechanism_comment"] += (
                            "Unsupported mechanism: " + record["mechanism"]
                        )
                    else:
                        valid_mechanism = record["mechanism"].lower().replace("-", " ")

                    final_draft["molecular_mechanism"] = {
                        "name": "",
                        "support": "inferred",
                    }
                    final_draft["mechanism_synopsis"] = []
                    # Mechanism evidence
                    final_draft["mechanism_evidence"] = []
                    # Save mechanism, mechanism evidence and comment in the source field
                    final_draft["source_data"]["mechanism"] = valid_mechanism
                    final_draft["source_data"]["mechanism_evidence"] = "; ".join(record["evidence"])
                    if "mechanism" in record["comments"]:
                        final_draft["source_data"]["mechanism_comment"] += "\n" + record["comments"]["mechanism"]

                    # Disease
                    disease_cross_references = []
                    if record["mondo_id"] != "":
                        disease_data = fetch_disease_info(api_url, record["mondo_id"])
                        disease_cross_references.append(
                            {
                                "source": "Mondo",
                                "identifier": record["mondo_id"],
                                "disease_name": disease_data["disease"].lower(),
                                "original_disease_name": disease_data["disease"],
                            }
                        )
                    if record["disease_id"] != "":
                        omim_list = record["disease_id"].split(",")
                        for omim_id in omim_list:
                            omim_id = omim_id.replace("OMIM:", "")
                            omim_id = omim_id.replace("MIM:", "")
                            if omim_id.strip().isdigit():
                                disease_data = fetch_disease_info(
                                    api_url, omim_id.strip()
                                )
                                disease_cross_references.append(
                                    {
                                        "source": "OMIM",
                                        "identifier": omim_id.strip(),
                                        "disease_name": disease_data["disease"].lower(),
                                        "original_disease_name": disease_data[
                                            "disease"
                                        ],
                                    }
                                )
                    final_draft["disease"] = {
                        "disease_name": "",
                        "cross_references": [],
                    }
                    # Save the disease name and cross references in the source field
                    final_draft["source_data"]["disease"] = record["disease"]
                    final_draft["source_data"]["disease_cross_references"] = disease_cross_references

                    final_draft["session_name"] = (
                        f"{source}_{record['gene_symbol']}_{final_draft['allelic_requirement']}"
                    )

                    # Check if the draft already exists based on:
                    # session name, gene, disease and allelic requirement
                    if final_draft["session_name"] in automatic_drafts:
                        existing_draft = automatic_drafts[final_draft["session_name"]]
                        if (
                            existing_draft["locus"] == final_draft["locus"]
                            and existing_draft["disease"]
                            == final_draft["disease"]["disease_name"]
                            and existing_draft["allelic_requirement"]
                            == final_draft["allelic_requirement"]
                        ):
                            print(
                                f"Draft record for session '{final_draft['session_name']}' already exists. Skipping insertion."
                            )
                            continue
                        else:
                            print(
                                f"Draft record for session '{final_draft['session_name']}' already exists but with different data. Consider reviewing it before inserting a new draft."
                            )
                            continue

                    # Call G2P API to insert the curation draft
                    insert_draft_record(
                         api_url,
                         cookies,
                         {"json_data": final_draft, "status": "automatic"},
                    )


def main():
    parser = argparse.ArgumentParser(description="Create draft records")
    parser.add_argument("--source", required=True, help="ClinGen or PanelApp")
    parser.add_argument(
        "--input_file",
        required=True,
        help="Data to be used to create drafts (supported formats: json and tsv)",
    )
    parser.add_argument(
        "--config", required=True, help="Config file with details to G2P API"
    )
    args = parser.parse_args()

    input_file = args.input_file
    source = args.source
    output_file = f"{source}_draft_records_created.json"

    # Load the config file
    config = configparser.ConfigParser()
    config.read(args.config)

    try:
        api = config["api"]
    except KeyError:
        sys.exit("ERROR: 'api' missing from config file")
    else:
        api_url = api["api_url"]
        api_username = api["api_username"]
        api_password = api["api_password"]

    if not os.path.isfile(input_file):
        sys.exit(f"Invalid input file '{input_file}'")

    cookies = login(api_username, api_password, api_url)

    # Get all existing automatic drafts to check for duplicates before inserting new ones
    automatic_drafts = get_all_automatic_drafts(api_url, cookies)

    prepare_draft_records(
        input_file, source.lower(), automatic_drafts, api_url, cookies
    )

    logout(api_url, cookies)


if __name__ == "__main__":
    main()
