#!/usr/bin/env python3

import os.path
import sys
import argparse
import re
import requests
import configparser


"""
    Script to merge G2P records.

    Params:
            --config : Config file name containing the API connection info (mandatory)
                    File format is the following:
                        [api]
                        api_url = <>

            --file : Tab delimited file with records to be merged (mandatory)
                File format is the following:
                    g2p id to keep\tg2p ids to merge\tdisease name\tgenotype

            --api_username: Username to connect to the G2P API (mandatory)
            --api_password: Password to connect to the G2P API (mandatory)
"""


def read_file(file: str, api_username: str, api_password: str, api_url: str) -> list:
    """
    Reads the records from the input file and creates the data to be sent to the API.
    The API runs most checks before merging the records.

    Args:
        file (str): tab delimited input file
        api_username (str):
        api_password (str):
        api_url (str):

    Returns:
        list: list of records to be merged

    Example output:
    [ 
        {"g2p_ids": ["G2P01201, G2P01203"], "final_g2p_id": "G2P00413"},
        {"g2p_ids": ["G2P01200"], "final_g2p_id": "G2P00001"}
    ]
    """
    file_not_merged = "records_not_merged.txt"
    records_to_merge = []

    with (open(file, "r") as fh, open(file_not_merged, "w") as wr):
        for line in fh:
            if not line.startswith("g2p id"):
                data = line.rstrip("\n").split("\t")
                # The data can be wrapped in quotes - remove them
                g2p_id_to_keep = data[0].strip().strip('"')
                g2p_ids_merge = [value.strip().strip('"') for value in data[1].split(",")]
                disease_name = data[2].strip()
                genotype = data[3].strip()

                if not re.search("^G2P\d{5}$", g2p_id_to_keep):
                    wr.write(line)

                else:
                    # Create output
                    records_to_merge.append({
                        "g2p_ids": g2p_ids_merge,
                        "final_g2p_id": g2p_id_to_keep
                    })

    return records_to_merge


def merge_records(records_to_merge: list, api_username: str, api_password: str, api_url: str) -> None:
    """
    Method to call the API to merge records.
    Endpoint: merge_records/

    Args:
        records_to_merge (list): list of records to merge
        api_username (str): G2P API username (super user)
        api_password (str): G2P API password
        api_url (str): URL to the G2P API
    """
    merge_url = "merge_records/"
    login_url = "login/"

    data = {"username": api_username, "password": api_password}

    response = requests.post(api_url + login_url, json=data)
    if response.status_code == 200:
        try:
            response_update = requests.post(
                api_url + merge_url, json=records_to_merge, cookies=response.cookies
            )
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
    else:
        print("Error: cannot login into G2P")


def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--config", required=True, help="Config file")
    parser.add_argument(
        "--file",
        required=True,
        help="Tab delimited file with all records to be merged",
    )
    parser.add_argument(
        "--api_username", required=True, help="Username to connect to the G2P API"
    )
    parser.add_argument(
        "--api_password", required=True, help="Password to connect to the G2P API"
    )
    args = parser.parse_args()

    file = args.file
    config_file = args.config
    api_username = args.api_username
    api_password = args.api_password

    # Load the config file
    config = configparser.ConfigParser()
    config.read(config_file)

    try:
        api_url = config["api"]["api_url"]
    except KeyError:
        sys.exit("Config: [api] is missing from the config file")

    if os.path.isfile(file):
        expected_columns = [
            "g2p id to keep",
            "g2p ids to merge",
            "disease name",
            "genotype"
        ]

        with open(file, "r", encoding="utf-8") as fh:
            header = fh.readline().strip().split("\t")
            if header[:4] != expected_columns:
                sys.exit(
                    f"Error: File format is incorrect. Found: {header[:4]}; Expected: {expected_columns}"
                )

        print("Parsing records to merge...")
        records_to_merge = read_file(file, api_username, api_password, api_url)
        print("Parsing records to merge... done\n")

        print("Merging records...")
        merge_records(records_to_merge, api_username, api_password, api_url)
        print("Merging records... done\n")

    else:
        print(f"Input file is invalid '{file}'")


if __name__ == "__main__":
    main()
