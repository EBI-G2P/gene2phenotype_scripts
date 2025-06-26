#!/usr/bin/env python3

import subprocess
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


def fetch_g2p_records(data: dict[str, Any]) -> requests.Response:
    """Fetches the G2P record using all the panels download
    User is not authenticated to ensure only visible panels is downloaded

    Args:
        data (dict[str, Any]): The DB and API configuration dictionary

    Returns:
        requests.Response: The response object from the get request
    """

    fetch_g2p_records = "panel/all/download/"

    url = data["api_url"] + fetch_g2p_records
    response = requests.get(url)

    if response.status_code == 200:
        return response

    else:
        sys.exit("Issues downloading the file")


def reading_data(response: requests.Response) -> list:
    """Reads the data from the response and converts it to a csv.DictReader

    Args:
        response (requests.Response): The response object from the get request

    Returns:
        list: The list reader
    """

    csv_content = io.StringIO(response.content.decode("utf-8"))
    reader = csv.DictReader(csv_content)

    return reader


def get_unsubmitted_record(data: dict[str, Any]) -> list:
    """Gets the unsubmitted records from the API

    Args:
        data (dict[str, Any]): DB and API configuration dictionary

    Returns:
        list: Returns a list of unsubmitted ids
    """
    fetch_unsubmitted_record = "unsubmitted_stable_ids/"

    url = data["api_url"] + fetch_unsubmitted_record
    response = requests.get(url)
    if response.status_code == 200:
        unsubmitted = json.loads(response.content)
        return unsubmitted
    else:
        print("Could not fetch unsubmitted records")


def get_later_review_date(data: dict[str, Any])-> list:
    """Gets later review date

    Args:
        data (dict[str, Any]): DB/api configuration

    Returns:
        list: Containing submission ids 
    """    
    fetch_later_review_date = "later_review_date/"

    url = data["api_url"] + fetch_later_review_date
    response = requests.get(url)
    if response.status_code == 200:
        later_date = json.loads(response.content)
        return later_date["ids"]
    else:
        print("Could not fetch later review date")
        

def retrieve_unsubmitted_records(data: dict[str, Any], unsubmitted: list) -> list:
    """Retrieves unsubmitted records from the read data from the file

    Args:
        data (dict[str, Any]): Content of the CSV file
        unsubmitted (list): List of unsubmitted ids

    Returns:
        list: List of the records that have not been submitted
    """
    records = [row for row in data if row["g2p id"] in unsubmitted]

    return records

def get_stable_id_associated_with_the_submission(data: dict[str, Any], submission_id: str) -> str:
    """Gets stable ids associated with the submissions

    Args:
        data (dict[str, Any]): DB/API configuration
        submission_id (str): Submission id

    Returns:
        str: Stable id strings which would be used in comparison
    """    

    fetch_record_data = f"submissions/{submission_id}"

    url = data["api_url"] + fetch_record_data

    response = requests.get(url)
    if response.status_code == 200:
        stable_id = json.loads(response.content)
        return stable_id[0]
    else:
        print("Could not get stable id associated with this submission")


def read_from_old_gencc_submission(file: str) -> list:
    """Read from old GenCC submission txt file

    Args:
        file (str): File name

    Returns:
        list: The list of data
    """    
    with open(file, "r", encoding="utf-8") as opened:
        reader = csv.DictReader(opened, delimiter='\t')
        data = list(reader)

    return data

def compare_data_changes(old_reader: list, later_date_ids: list, new_reader: list, data: dict[str, Any])-> list:
    """Compare data changes between what has been submitted and what is about to be unsubmitted for submission with later review date

    Args:
        old_reader (list): Old file data
        later_date_ids (list): The later date ids
        new_reader (list): New G2P file data
        data (dict[str, Any]): DB/API configuration

    Returns:
        list: List of data that has had changes.
    """    
    updated_data = []

    new_data_lookup = {row["g2p id"]: row for row in new_reader}
    for old_row in old_reader:
        submission_id = old_row["submission_id"]
        if submission_id not in later_date_ids:
            continue

        stable_id = get_stable_id_associated_with_the_submission(data, submission_id)
        new_row = new_data_lookup.get(stable_id)

        if not new_row:
            continue  
        
        if (
            new_row.get("disease mim") != old_row.get("disease id")
            and new_row.get("disease MONDO") != old_row.get("disease id")
        ):  
            new_row["update"] = 1
            updated_data.append(new_row)
        elif new_row.get("disease name") != old_row.get("disease name"):
            new_row["update"] = 1
            updated_data.append(new_row)
        elif new_row.get("allelic requirement") != old_row.get("moi_name"):
            new_row["update"] = 1
            updated_data.append(new_row)
        elif new_row.get("confidence category") != old_row.get("classification_name"):
            new_row["update"] = 1
            updated_data.append(new_row)
        elif new_row.get("publications") != old_row.get("pmids"):
            new_row["update"] = 1
            updated_data.append(new_row)

    return updated_data

def create_gencc_submission_record( submission_id: str, date: str, type_of: str, stable_id: str
):
    """Creates the GenCC submission record by creating the record in the gencc_submission table of the DB

    Args:
        submission_id (str): Submission id
        date (str): date
        type_of (str): Usually create
        stable_id (str): G2P stable id
    """
    create_info = {
        "submission_id": submission_id,
        "date_of_submission": date,
        "type_of_submission": type_of,
        "g2p_stable_id": stable_id,
    }

    return create_info

def post_gencc_submission(list_of_data: list, data:dict[str, Any]):
    """Create GenCC submission record, using bulk create

    Args:
        list_of_data (list): _description_
        data (dict[str, Any]): _description_

    """
    
    create_url = "gencc_create/"

    login_url = "login/"

    login_info = {"username": data["username"], "password": data["api_password"]}

    response = requests.post(data["api_url"] + login_url, json=login_info)
    if response.status_code == 200:
        try:
            response_create = requests.post(
                data["api_url"] + create_url, json=list_of_data, cookies=response.cookies
            )
            if response_create.status_code in (200, 201):
                print(
                    f"GenCC submission for the records was created successfully"
                )
            else:
                print(
                    f"Issues creating the submissions {response_create.status_code} {response_create.json}"
                )
        except Exception as e:
            print("Error:", e)

def add_unsubmitted_ids_and_later_review_date(updated_data: list, data: list)-> list:
    """Merging data

    Args:
        updated_data (list): Updated data from later review date checks
        data (list): Unsubmitted ids data

    Returns:
        list: Merged list of both of them
    """    
    return updated_data + data

def write_to_the_GenCC_file(
    data: list, outfile: str, dry: str, db_config: dict[str, Any], type_of: str=None
) -> str:
    """Creates the G2P_GenCC.txt and also calls the create the gencc_submission function when dry is False,
    A real run
    Also writes out records with no disease id

    Args:
        data (dict[str, Any]): Record of unsubmitted ids
        outfile (str): Txt file to be created
        dry (str): To allow for a real run, which updates the db
        db_config (dict[str, Any]): DB/API configuration dictionary

    Returns:
        str: An output file
    """
    with open(outfile, mode="w") as output_file:
        output_file.write(
            "submission_id\thgnc_id\thgnc_symbol\tdisease_id\tdisease_name\tmoi_id\tmoi_name\tsubmitter_id\tsubmitter_name\tclassification_id\tclassification_name\tdate\tpublic_report_url\tpmids\tassertion_criteria_url\n"
        )
        issues_with_record = []
        gencc_list = []
        submission_id_base = "1000112"
        for record in data:
            g2p_id = record["g2p id"]
            submission_id = submission_id_base + str(g2p_id[3:])

            hgnc_id = record["hgnc id"]
            hgnc_symbol = record["gene symbol"]

            disease_id = record["disease mim"] or record["disease MONDO"]
            if disease_id is None or disease_id == "":
                issues_with_record.append(g2p_id)
                continue
            disease_name = record["disease name"]
            moi_id = allelic_requirement[record["allelic requirement"]]
            moi_name = record["allelic requirement"]
            submitter_id = "GENCC:000113"
            submitter_name = "G2P"
            classification_id = confidence_category[record["confidence"]]
            classification_name = record["confidence"]
            dt = datetime.fromisoformat(record["date of last review"])
            date = dt.strftime("%Y/%m/%d")
            record_url = "https://www.ebi.ac.uk/gene2phenotype/lgd/" + g2p_id
            pmids = record["publications"]
            assertion_criteria_url = "https://www.ebi.ac.uk/gene2phenotype/terminology"

            line_to_output = f"{submission_id}\t{hgnc_id}\t{hgnc_symbol}\t{disease_id}\t{disease_name}\t{moi_id}\t{moi_name}\t{submitter_id}\t{submitter_name}\t{classification_id}\t{classification_name}\t{date}\t{record_url}\t{pmids}\t{assertion_criteria_url}\n"
            output_file.write(line_to_output)
            db_date = create_datetime_now()
            if dry is True:
                if type_of:
                    type_of = type_of
                    created_record = create_gencc_submission_record(
                        submission_id, db_date, type_of, g2p_id
                    )
                    gencc_list.append(created_record)
                else:
                    update = record["update"]
                    if update:
                        type_of = "update"
                        created_record = create_gencc_submission_record(submission_id, db_date, type_of, g2p_id)
                        gencc_list.append(created_record)
                    else:
                        type_of = "create"
                        created_record = create_gencc_submission_record(submission_id, db_date, type_of, g2p_id)
                        gencc_list.append(created_record)
    if len(issues_with_record) > 0:
        with open("record_with_issues.txt", mode="w") as textfile:
            for issues in issues_with_record:
                textfile.write(f"{issues}\n")
    
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
    """Reads from Configuration file (config.ini) and creates a dictionary of data

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

def get_output_paths(path: str) -> Tuple[ str, str, str]:
    """Gets the output path from the args.path

    Args:
        path (str): The path (args.path)

    Returns:
        Tuple[str, str, str]: Gencc dir (args.path if given), output file and final_output_file
    """    
    if path:
        timestamp = create_datetime_now()
        gencc_dir = os.path.join(path, timestamp)
        os.makedirs(gencc_dir, exist_ok=True) # to make the directory

        output_file = os.path.join(gencc_dir, "G2P_GenCC.txt")
        final_output_file = os.path.join(gencc_dir, "G2P_GenCC.xlsx")
        return gencc_dir, output_file, final_output_file
    else:
        output_file = "G2P_GenCC.txt"
        final_output_file = "G2P_GenCC.xlsx"
        return ".", output_file, final_output_file

def handle_new_submission(read_data: list, output_file: str, dry: bool, db_config: dict[str, Any]) -> write_to_the_GenCC_file:
    """Handles new submission, so new generation of the G2P records and GenCC data 

    Args:
        read_data (list): Data from the G2P all records file
        output_file (str): Output text file where the submission will be written into
        dry (bool): Boolean (args.dry)
        db_config (dict[str, Any]): DB/API configuration

    Returns:
        write_to_the_GenCC_file: A text file containing all the records to be submitted to GenCC
    """
    type_of = "create"    
    return write_to_the_GenCC_file(read_data, output_file, dry, db_config, type_of)

def handle_existing_submission(read_data: list, output_file: str, dry: bool, db_config: dict[str, Any], old_file: str) -> write_to_the_GenCC_file:
    """Handles existing submission.

    Args:
        read_data (list): Data from the G2P all records file
        output_file (str): Output text file where the submission will be written into
        dry (bool): Boolean (args.dry)
        db_config (dict[str, Any]): DB/API configuration
        old_file (str) : Old file containing submitted GenCC submission to allow comparison

    Returns:
        write_to_the_GenCC_file: A text file containing all the records to be submitted to GenCC
    """    
    unsubmitted = get_unsubmitted_record(db_config)
    common = retrieve_unsubmitted_records(read_data, unsubmitted)

    if old_file:
        old_reader = read_from_old_gencc_submission(old_file)
        later_review = get_later_review_date(db_config)
        compared = compare_data_changes(old_reader, later_review, read_data, db_config)
        merged_data = add_unsubmitted_ids_and_later_review_date(compared, common)
        return write_to_the_GenCC_file(merged_data, output_file, dry, db_config)
    else:
        type_of = "create"
        return write_to_the_GenCC_file(common, output_file, dry, db_config, type_of)

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
        "--write_to_db",
        required=False,
        action="store_true",
    )
    ap.add_argument(
        "--new",
        action="store_true",
        required=False,
        help="New is used to determine if this is a new or clean GenCC submission, so no checks on the submission id or stable id will be done",

    )
    ap.add_argument("--old_file", required=False, help="Old file to compare changes in G2P ids in txt format")
    args = ap.parse_args()

    if args.new and args.old_file:
        ap.error("Cannot use --new and --old_file at the same time.")

    db_config = read_from_config_file(args.config_file)
    dry_run = args.write_to_db

    gencc_dir, output_file, final_output_file = get_output_paths(args.path)

    print("Downloading G2P files...")
    file_data = fetch_g2p_records(db_config)
    print("Reading data from downloaded file")
    read_data = reading_data(file_data)

    if args.new:
        print("Handling new submission")
        outfile, gencc_list = handle_new_submission(read_data, output_file, dry_run, db_config)
    else:
        print("Handling existing submission")
        outfile, gencc_list = handle_existing_submission(read_data, output_file, dry_run, db_config, args.old_file)

    print("Converting text file to Excel file")
    convert_txt_to_excel(outfile, final_output_file)

    if args.write_to_db:
        print("Writing submission to the database")
        post_gencc_submission(gencc_list, db_config)


if __name__ == "__main__":
    main()
