#!/usr/bin/env python3

import argparse
import os
import sys
from datetime import date, datetime
import requests
import json
import csv
from openpyxl import Workbook
from typing import Any, Tuple
from configparser import ConfigParser
import io

# Mapping terms to GenCC IDs
allelic_requirement = {
    "biallelic_autosomal": "HP:0000007",
    "monoallelic_autosomal": "HP:0000006",
    "monoallelic_X_hemizygous": "HP:0001417",
    "monoallelic_Y_hemizygous": "HP:0001450",
    "monoallelic_X_heterozygous": "HP:0001417",
    "mitochondrial": "HP:0001427",
    "monoallelic_PAR": "HP:0000006",
    "biallelic_PAR": "HP:0000007",
}
confidence_category = {
    "definitive": "GENCC:100001",
    "strong": "GENCC:100002",
    "moderate": "GENCC:100003",
    "limited": "GENCC:100004",
    "disputed": "GENCC:100005",
    "refuted": "GENCC:100006",
}


def fetch_g2p_records(data: dict[str, Any]) -> list:
    """
    Fetches the G2P record using download for all panels.
    User is not authenticated to ensure only visible panels are downloaded

    Args:
        data (dict[str, Any]): The DB and API configuration dictionary

    Returns:
        list: The list reader
    """
    url = f"{data['api_url'].rstrip('/')}/panel/all/download/"
    response = requests.get(url)

    if response.status_code == 200:
        csv_content = io.StringIO(response.content.decode("utf-8"))
        reader = csv.DictReader(csv_content)
        return list(reader)

    else:
        sys.exit(f"Issues downloading the G2P file from {url}")


def get_unsubmitted_records(data: dict[str, Any]) -> list:
    """
    Gets the G2P records that haven't been submitted to GenCC yet.

    Args:
        data (dict[str, Any]): DB and API configuration dictionary

    Returns:
        list: Returns a list of unsubmitted record ids
    """
    url = f"{data['api_url'].rstrip('/')}/unsubmitted_stable_ids/"
    response = requests.get(url)
    if response.status_code == 200:
        unsubmitted = json.loads(response.content)
        return unsubmitted
    else:
        sys.exit(f"Could not fetch unsubmitted records from {url}")


def get_updated_records(data: dict[str, Any]) -> dict:
    """
    Gets G2P records that have been udpated since last GenCC submission

    Args:
        data (dict[str, Any]): DB/API configuration

    Returns:
        dict: Dictionary containing records that have been udpated since last submission
    """
    url = f"{data['api_url'].rstrip('/')}/later_review_date/"
    response = requests.get(url)
    if response.status_code == 200:
        response_data = json.loads(response.content)
        return response_data
    else:
        sys.exit(f"Could not fetch updated records from {url}")


def get_deleted_records(data: dict[str, Any]) -> dict:
    """
    Gets G2P records that have been submitted to GenCC but are now deleted in G2P

    Args:
        data (dict[str, Any]): DB/API configuration

    Returns:
        dict: Dictionary containing records that have been udpated since last submission
    """
    url = f"{data['api_url'].rstrip('/')}/gencc_deleted_records/"
    response = requests.get(url)
    if response.status_code == 200:
        response_data = json.loads(response.content)
        return response_data
    else:
        sys.exit(f"Could not fetch deleted records from {url}")


def post_gencc_submission(list_of_data: list, db_config: dict[str, Any]):
    """
    Calls the G2P API to write to the db which records are going to be submitted to GenCC.

    Args:
        list_of_data (list): list of records to write to the db
        db_config (dict[str, Any]): dictionary with the database connection details and the API URL

    """
    create_url = "gencc_create/"
    login_url = "login/"

    login_info = {
        "username": db_config["username"],
        "password": db_config["api_password"],
    }

    response = requests.post(db_config["api_url"] + login_url, json=login_info)
    if response.status_code == 200:
        try:
            response_create = requests.post(
                db_config["api_url"] + create_url,
                json=list_of_data,
                cookies=response.cookies,
            )
            if response_create.status_code in (200, 201):
                print(f"GenCC submission for the records was created successfully")
            else:
                sys.exit(
                    f"Issues creating the submissions {response_create.status_code} {response_create.json}"
                )
        except Exception as e:
            sys.exit("Error:", e)


def write_to_the_GenCC_file(
    records_data: list,
    outfile: str
) -> str:
    """
    Write the records to be submitted to the output file 'G2P_GenCC.txt'.
    Records without disease id are excluded from the submission and written
    to file 'record_with_issues.txt'.

    Args:
        records_data (dict[str, Any]): Record of unsubmitted ids
        outfile (str): Txt file to be created

    Returns:
        outfile (str): The output file with records to be submitted
        gencc_list (list): list of records to write to the db
    """
    with open(outfile, mode="w") as rw:
        rw.write(
            "submission_id\thgnc_id\thgnc_symbol\tdisease_id\tdisease_name\tmoi_id\tmoi_name"
            "\tsubmitter_id\tsubmitter_name\tclassification_id\tclassification_name\tdate"
            "\tpublic_report_url\tnotes\tpmids\tassertion_criteria_url\n"
        )
        issues_with_record = {}
        gencc_list = []
        submission_id_base = "1000112"
        submitter_id = "GENCC:000112"
        submitter_name = "TGMI"  # To be updated in the future
        assertion_criteria_url = "https://www.ebi.ac.uk/gene2phenotype/about/terminology"

        for record in records_data:
            g2p_id = record["g2p id"]

            # Prepare submission ID
            # Create a new submission ID if it's a new record
            if "type_of_submission" not in record or record["type_of_submission"] == "create":
                submission_id = submission_id_base + str(g2p_id[3:]) # TODO: review
                type_of_submission = "create" # new records don't always have the key 'type_of_submission'
            else:
                submission_id = record["submission_id"]
                type_of_submission = record["type_of_submission"]

            hgnc_id = record["hgnc id"]
            hgnc_symbol = record["gene symbol"]

            disease_id = record["disease mim"] or record["disease MONDO"]
            # We cannot submit records without disease ID
            if disease_id is None or disease_id == "":
                issues_with_record[g2p_id] = "Missing disease ID"
                continue

            disease_name = record["disease name"]

            if record["allelic requirement"] in allelic_requirement:
                moi_id = allelic_requirement[record["allelic requirement"]]
            else:
                issues_with_record[g2p_id] = f"Unsupported allelic requirement '{record['allelic requirement']}'"
                continue

            moi_name = record["allelic requirement"]
            classification_id = confidence_category[record["confidence"]]
            classification_name = record["confidence"]
            dt = datetime.fromisoformat(record["date of last review"])
            date = dt.strftime("%Y/%m/%d")
            record_url = "https://www.ebi.ac.uk/gene2phenotype/lgd/" + g2p_id
            pmids = record["publications"]

            line_to_output = f"{submission_id}\t{hgnc_id}\t{hgnc_symbol}\t{disease_id}\t{disease_name}\t{moi_id}\t{moi_name}\t{submitter_id}\t{submitter_name}\t{classification_id}\t{classification_name}\t{date}\t{record_url}\t\t{pmids}\t{assertion_criteria_url}\n"
            rw.write(line_to_output)
            db_date = create_datetime_now()
            
            gencc_list.append(
                {
                    "submission_id": submission_id,
                    "date_of_submission": db_date,
                    "type_of_submission": type_of_submission,
                    "g2p_stable_id": g2p_id,
                }
            )

    if issues_with_record:
        with open("record_with_issues.txt", mode="w") as textfile:
            for g2p_id in issues_with_record:
                textfile.write(f"{g2p_id}\t{issues_with_record[g2p_id]}\n")

    return outfile, gencc_list


def create_datetime_now() -> date:
    """Creates datetime at present and formats it to YYYY-MM-DD

    Returns:
        date: Formatted date
    """

    time_now = datetime.now()
    db_date = time_now.strftime("%Y-%m-%d")
    return db_date


def convert_txt_to_excel(input_file: str, output_file: str):
    """Converts txt to Excel file

    Args:
        input_file (str): Text file to be converted
        output_file (str): Excel file that will be created
    """
    wb = Workbook()
    ws = wb.active

    delimiter = "\t"

    with open(input_file, "r") as file:
        for row_index, line in enumerate(file, start=1):
            for col_index, value in enumerate(line.strip().split(delimiter), start=1):
                ws.cell(row=row_index, column=col_index, value=value)
    wb.save(output_file)


def read_from_config_file(config_file: str) -> dict[str, Any]:
    """
    Reads the database and api info from the config file.

    Args:
        config_file (str): configuration file

    Returns:
        dict[str, Any]: DB/API configuration dict
    """
    data = {}
    config = ConfigParser()
    config.read(config_file)
    data["host"] = config["database"]["host"]
    data["port"] = config["database"]["port"]
    data["user"] = config["database"]["user"]
    data["database"] = config["database"]["name"]
    data["password"] = config["database"]["password"]
    data["api_url"] = config["api"]["api_url"]
    data["username"] = config["api"]["username"]
    data["api_password"] = config["api"]["password"]

    return data


def get_output_paths(path: str) -> Tuple[str, str]:
    """
    Create the output directory and prepare the output file path.

    Args:
        path (str): path where the G2P and GenCC files are going to be saved

    Returns:
        Tuple[str, str]: output file (txt format) and final_output_file (xlsx format)
    """
    if path:
        timestamp = create_datetime_now()
        gencc_dir = os.path.join(path, timestamp)
        os.makedirs(gencc_dir, exist_ok=True)  # to make the directory

        output_file = os.path.join(gencc_dir, "G2P_GenCC.txt")
        final_output_file = os.path.join(gencc_dir, "G2P_GenCC.xlsx")
    else:
        # TODO: this should create the new directory too
        output_file = "G2P_GenCC.txt"
        final_output_file = "G2P_GenCC.xlsx"

    return output_file, final_output_file


def handle_existing_submission(
    g2p_data: list,
    output_file: str,
    db_config: dict[str, Any],
) -> None:
    """Handles existing submission.

    Args:
        g2p_data (list): Data from the G2P all records file
        output_file (str): Output text file where the submission will be written into
        db_config (dict[str, Any]): DB/API configuration
    """
    # Retrieves the records that have been updated since last submission
    # We only need to re-submit these if the disease name/IDs have changed
    updated_records_data = get_updated_records(db_config)
    # TODO: review the updated records

    # Retrieves the records that have been deleted in G2P but were submitted to GenCC
    deleted_records_data = get_deleted_records(db_config)

    # List of records that haven't been submitted yet
    # These records are part of public panels and are live in G2P (not deleted)
    unsubmitted = get_unsubmitted_records(db_config)
    # Retrieves the unsubmitted records data from the G2P downloaded file
    records_to_submit = []

    for row in g2p_data:
        # Prepare unsubmitted records
        if row["g2p id"] in unsubmitted:
            row["type_of_submission"] = "create"
            row["submission_id"] = None
            records_to_submit.append(row)
        if row["g2p id"] in list(updated_records_data["ids"].keys()):
            row["type_of_submission"] = "update"
            row["submission_id"] = updated_records_data["ids"][row["g2p id"]]
            records_to_submit.append(row)

    outfile, gencc_list = write_to_the_GenCC_file(records_to_submit, output_file)

    # TODO: write deleted records to file

    return outfile, gencc_list


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "-p",
        "--path",
        help="Path where the G2P and GenCC files are going to be saved",
        required=False,
    )
    ap.add_argument("--config_file", required=True, help="G2P Configuration file")
    ap.add_argument(
        "--skip_write_to_db",
        required=False,
        action="store_true",
    )
    ap.add_argument(
        "--new_submission",
        action="store_true",
        required=False,
        help="New is used to determine this is a new submission to GenCC",
    )
    args = ap.parse_args()

    db_config = read_from_config_file(args.config_file)
    skip_write_to_db = args.skip_write_to_db

    output_file, final_output_file = get_output_paths(args.path)

    print("Getting G2P records from the API")
    g2p_data = fetch_g2p_records(db_config)

    if args.new_submission:
        print("Handling new submission")
        outfile, gencc_list = write_to_the_GenCC_file(g2p_data, output_file)
    else:
        print("Handling existing submission")
        outfile, gencc_list = handle_existing_submission(g2p_data, output_file, db_config)

    # print("Converting text file to Excel file")
    # convert_txt_to_excel(outfile, final_output_file)

    # if not skip_write_to_db:
    #     print("Writing submission details to the database")
    #     post_gencc_submission(gencc_list, db_config)


if __name__ == "__main__":
    main()
