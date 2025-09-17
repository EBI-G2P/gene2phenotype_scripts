#!/usr/bin/env python3

import csv
import sys
import argparse
import requests
import configparser
from pathlib import Path
from typing import Dict, List


"""
    Script to merge G2P records.

    Params:
            --config : Config file name containing the API connection info (mandatory)
                    File format is the following:
                        [api]
                        api_url = <>

            --input-file : Tab delimited file with records to be merged (mandatory)
                File format is the following:
                    g2p id to keep\tg2p ids to merge\tgene\tdisease name\tgenotype\tmechanism

            --api-username: Username to connect to the G2P API (mandatory)
            --api-password: Password to connect to the G2P API (mandatory)
            --dryrun: Test script without running the updates (default: 0)
"""


EXPECTED_COLUMNS = [
    "g2p id to keep",
    "g2p ids to merge",
    "gene",
    "disease name",
    "genotype",
    "mechanism"
]


def load_records(file: Path) -> List[Dict[str, str]]:
    """Read the input text file and validate the header."""
    with file.open(newline="", encoding='utf-8') as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        header = reader.fieldnames or []
        if header[: len(EXPECTED_COLUMNS)] != EXPECTED_COLUMNS:
            raise ValueError(
                f"File format is incorrect. Found: {header[:len(EXPECTED_COLUMNS)]}; "
                f"expected: {EXPECTED_COLUMNS}"
            )
        return list(reader)


def login(api_username: str, api_password: str, api_url: str) -> requests.cookies.RequestsCookieJar:
    """Login into G2P API"""
    login_url = f"{api_url.rstrip('/')}/login/"

    response = requests.post(login_url, json={"username": api_username, "password": api_password})

    if response.status_code != 200:
        sys.exit("Login failed. Check your credentials and API URL.")

    return response.cookies


def logout(api_url: str, cookies: requests.cookies.RequestsCookieJar) -> None:
    """Logout of the API"""
    logout_url = f"{api_url.rstrip('/')}/logout/"

    response = requests.post(logout_url, cookies=cookies)

    if response.status_code != 204:
        sys.exit("Logout failed. Check your credentials and API URL.")


def process_records(records: List[Dict[str, str]], api_url: str, cookies: requests.cookies.RequestsCookieJar) -> list:
    """
    Process the records and creates the data to be sent to the API.
    It checks if the gene, disease, genotype and mechanism from the file match the record to keep.
    All the other checks before merging the records are done by the API.

    Args:
        records (List[Dict[str, str]]): records from the input file to be processed
        api_url (str): URL to the G2P API
        cookies (requests.cookies.RequestsCookieJar): cookies contained tokens (login credentials)

    Returns:
        list: list of records to be merged

    Example output:
    [ 
        {"g2p_ids": ["G2P01201, G2P01203"], "final_g2p_id": "G2P00413"},
        {"g2p_ids": ["G2P01200"], "final_g2p_id": "G2P00001"}
    ]
    """
    file_not_merged = "records_not_merged.txt"
    file_to_merge = "records_to_merge.txt"
    records_to_merge = []

    lgd_url = f"{api_url.rstrip('/')}/lgd/"

    with (open(file_not_merged, "w") as wr, open(file_to_merge, "w") as wr_merge):
        for line in records:
            g2p_id_to_keep = line["g2p id to keep"].strip()
            g2p_ids_merge = [x.strip() for x in line["g2p ids to merge"].split(',')]
            gene = line["gene"].strip()
            disease_name = line["disease name"].strip()
            genotype = line["genotype"].strip()
            mechanism = line["mechanism"].strip()

            if not g2p_id_to_keep.startswith("G2P"):
                wr.write(f"Invalid G2P ID\t{g2p_id_to_keep}\t{line['g2p ids to merge']}\n")
                continue

            try:
                lgd_response = requests.get(lgd_url+g2p_id_to_keep, cookies=cookies)
            except Exception as e:
                print(f"Error while fetching record {g2p_id_to_keep}:", str(e))
            else:
                if lgd_response.status_code == 200:
                    response_json = lgd_response.json()
                    if gene != response_json["locus"]["gene_symbol"]:
                        wr.write(f"Gene doesn't match\t{g2p_id_to_keep}\t{line['g2p ids to merge']}\n")
                    elif disease_name != response_json["disease"]["name"]:
                        wr.write(f"Disease doesn't match\t{g2p_id_to_keep}\t{line['g2p ids to merge']}\n")
                    elif genotype != response_json["genotype"]:
                        wr.write(f"Genotype doesn't match\t{g2p_id_to_keep}\t{line['g2p ids to merge']}\n")
                    elif mechanism != response_json["molecular_mechanism"]["mechanism"]:
                        wr.write(f"Mechanism doesn't match\t{g2p_id_to_keep}\t{line['g2p ids to merge']}\n")
                    else:
                        # Create output
                        records_to_merge.append({
                            "g2p_ids": g2p_ids_merge,
                            "final_g2p_id": g2p_id_to_keep
                        })
                        wr_merge.write(f"Merge {line['g2p ids to merge']} into {g2p_id_to_keep}\n")
                else:
                    print(f"Failed fetching record {g2p_id_to_keep}. Status code: {lgd_response.status_code}")

    return records_to_merge


def merge_records(records_to_merge: list, api_url: str, cookies:requests.cookies.RequestsCookieJar) -> None:
    """
    Method to call the API to merge records.
    Endpoint: merge_records/

    Args:
        records_to_merge (list): list of records to merge
        api_url (str): URL to the G2P API
        cookies (requests.cookies.RequestsCookieJar): cookies contained tokens (login credentials)
    """
    merge_url = f"{api_url.rstrip('/')}/merge_records/"

    try:
        response_update = requests.post(merge_url, json=records_to_merge, cookies=cookies)
    except Exception as e:
        print("Error:", e)
    else:
        if response_update.status_code == 200:
            response_json = response_update.json()
            print("Records merged successfully:", response_json)
            if "error" in response_json:
                print(
                "Failed to merge the following records:",
                response_json.get("error"),
            )
        elif response_update.status_code == 400 and "error" in response_update.json():
            print(
                "Failed to merge the following records:",
                response_update.json().get("error"),
            )
        else:
            print(
                f"Failed to merge records. Status code: {response_update.status_code}. Message:",
                response_update.json(),
            )


def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--config", required=True, help="Config file")
    parser.add_argument(
        "--input-file",
        required=True,
        help="Tab delimited file with all records to be merged",
    )
    parser.add_argument(
        "--api-username", required=True, help="Username to connect to the G2P API"
    )
    parser.add_argument(
        "--api-password", required=True, help="Password to connect to the G2P API"
    )
    parser.add_argument(
        "--dryrun", action="store_true", help="Option to test the which records are going to be merged"
    )
    args = parser.parse_args()

    input_file = args.input_file
    config_file = args.config
    api_username = args.api_username
    api_password = args.api_password
    dryrun = args.dryrun

    # Load the config file
    config = configparser.ConfigParser()
    config.read(config_file)

    try:
        api_url = config["api"]["api_url"]
    except KeyError:
        sys.exit("Config: [api] is missing from the config file")

    input_path = Path(input_file)
    if not input_path.is_file():
        sys.exit(f"Input file not found: {input_path}")

    try:
        # Load the records from the input file
        records = load_records(input_path)
    except ValueError as e:
        sys.exit(e)

    print("Logging in...")
    cookies = login(api_username, api_password, api_url)

    print("Parsing records to merge...")
    records_to_merge = process_records(records, api_url, cookies)
    print("Parsing records to merge... done")

    if not dryrun:
        print("Merging records...")
        merge_records(records_to_merge, api_url, cookies)
        print("Merging records... done")
    else:
        print("Dry run: merge skipped")

    print("Logging out...")
    logout(api_url, cookies)


if __name__ == "__main__":
    main()
