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


def fetch_g2p_records(data: dict[str, Any]) -> requests.Response:
    """
    Fetches the G2P record using all the panels download
    User is not authenticated to ensure only visible panels are downloaded

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
    """
    Reads the data from the response and converts it to a csv.DictReader

    Args:
        response (requests.Response): The response object from the get request

    Returns:
        list: The list reader
    """
    csv_content = io.StringIO(response.content.decode("utf-8"))
    reader = csv.DictReader(csv_content)

    return list(reader)


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
        sys.exit("Could not fetch unsubmitted records from API")


def get_later_review_date(data: dict[str, Any]) -> list:
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
        sys.exit("Could not fetch later review date from API")


def get_stable_id_associated_with_the_submission(
    data: dict[str, Any], submission_id: str
) -> str:
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
        sys.exit(
            f"Could not get the G2P stable ID associated with submission {submission_id}"
        )


def read_from_old_gencc_submission(file: str) -> list:
    """Read from old GenCC submission txt file

    Args:
        file (str): File name

    Returns:
        list: The list of data
    """
    with open(file, "r", encoding="utf-8") as opened:
        reader = csv.DictReader(opened, delimiter="\t")
        data = list(reader)

    return data


def compare_data_changes(
    old_reader: list, later_date_ids: list, new_reader: list, data: dict[str, Any]
) -> list:
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

        if new_row.get("disease mim") != old_row.get("disease id") and new_row.get(
            "disease MONDO"
        ) != old_row.get("disease id"):
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
    data: list, outfile: str, write_to_db: bool, type_of: str = None
) -> str:
    """
    Write the records to be submitted to the output file 'G2P_GenCC.txt'.
    Records without disease id are excluded from the submission and written
    to file 'record_with_issues.txt'.

    Args:
        data (dict[str, Any]): Record of unsubmitted ids
        outfile (str): Txt file to be created
        write_to_db (bool): Writes submission details to the db (table 'gencc_submission')
        type_of (str): If not provided, it is set to None

    Returns:
        outfile (str): The output file with records to be submitted
        gencc_list (list): list of records to write to the db
    """
    with open(outfile, mode="w") as output_file:
        output_file.write(
            "submission_id\thgnc_id\thgnc_symbol\tdisease_id\tdisease_name\tmoi_id\tmoi_name\tsubmitter_id\tsubmitter_name\tclassification_id\tclassification_name\tdate\tpublic_report_url\tnotes\tpmids\tassertion_criteria_url\n"
        )
        issues_with_record = []
        gencc_list = []
        submission_id_base = "1000112"
        for record in data:
            g2p_id = record["g2p id"]
            submission_id = (
                submission_id_base + str(g2p_id[3:])
                if not record.get("submission_id")
                else record.get("submission_id")
            )

            hgnc_id = record["hgnc id"]
            hgnc_symbol = record["gene symbol"]

            disease_id = record["disease mim"] or record["disease MONDO"]
            if disease_id is None or disease_id == "":
                issues_with_record.append(g2p_id)
                continue
            disease_name = record["disease name"]
            moi_id = allelic_requirement[record["allelic requirement"]]
            moi_name = record["allelic requirement"]
            submitter_id = "GENCC:000112"
            submitter_name = "TGMI"  # To be updated in the future
            classification_id = confidence_category[record["confidence"]]
            classification_name = record["confidence"]
            dt = datetime.fromisoformat(record["date of last review"])
            date = dt.strftime("%Y/%m/%d")
            record_url = "https://www.ebi.ac.uk/gene2phenotype/lgd/" + g2p_id
            pmids = record["publications"]
            assertion_criteria_url = (
                "https://www.ebi.ac.uk/gene2phenotype/about/terminology"
            )

            line_to_output = f"{submission_id}\t{hgnc_id}\t{hgnc_symbol}\t{disease_id}\t{disease_name}\t{moi_id}\t{moi_name}\t{submitter_id}\t{submitter_name}\t{classification_id}\t{classification_name}\t{date}\t{record_url}\t\t{pmids}\t{assertion_criteria_url}\n"
            output_file.write(line_to_output)
            db_date = create_datetime_now()
            if write_to_db:
                if not type_of:
                    # TODO: fix bug
                    type_of = "update" if record.get("update") == 1 else "create"
                gencc_list.append(
                    {
                        "submission_id": submission_id,
                        "date_of_submission": db_date,
                        "type_of_submission": type_of,
                        "g2p_stable_id": g2p_id,
                    }
                )
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
    read_data: list,
    output_file: str,
    write_to_db: bool,
    db_config: dict[str, Any],
    old_file: str,
) -> write_to_the_GenCC_file:
    """Handles existing submission.

    Args:
        read_data (list): Data from the G2P all records file
        output_file (str): Output text file where the submission will be written into
        write_to_db (bool): If set to True, write submission details to the db
        db_config (dict[str, Any]): DB/API configuration
        old_file (str) : Old file containing submitted GenCC submission to allow comparison

    Returns:
        write_to_the_GenCC_file: A text file containing all the records to be submitted to GenCC
    """
    unsubmitted = get_unsubmitted_record(db_config)
    # Retrieves unsubmitted records from the data from the file
    common = [row for row in read_data if row["g2p id"] in unsubmitted]

    if old_file:
        old_reader = read_from_old_gencc_submission(old_file)
        later_review = get_later_review_date(db_config)
        compared = compare_data_changes(old_reader, later_review, read_data, db_config)
        merged_data = compared + common
        return write_to_the_GenCC_file(merged_data, output_file, write_to_db)
    else:
        return write_to_the_GenCC_file(common, output_file, write_to_db, "create")


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
    ap.add_argument(
        "--old_file",
        required=False,
        help="Old file to compare changes in G2P ids in txt format",
    )
    args = ap.parse_args()

    if args.new and args.old_file:
        ap.error("Cannot use --new and --old_file at the same time.")

    db_config = read_from_config_file(args.config_file)
    write_to_db = args.write_to_db
    old_file = args.old_file

    if old_file and not os.path.isfile(old_file):
        sys.exit(f"Invalid file '{old_file}'")

    output_file, final_output_file = get_output_paths(args.path)

    print("Downloading G2P files")
    file_data = fetch_g2p_records(db_config)

    print("Reading data from downloaded file")
    read_data = reading_data(file_data)

    if args.new:
        print("Handling new submission")
        outfile, gencc_list = write_to_the_GenCC_file(
            read_data, output_file, write_to_db, "create"
        )
    else:
        print("Handling existing submission")
        outfile, gencc_list = handle_existing_submission(
            read_data, output_file, write_to_db, db_config, old_file
        )

    print("Converting text file to Excel file")
    convert_txt_to_excel(outfile, final_output_file)

    if write_to_db:
        print("Writing submission details to the database")
        post_gencc_submission(gencc_list, db_config)


if __name__ == "__main__":
    main()
