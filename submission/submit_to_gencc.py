#!/usr/bin/env python3

import argparse
import os
import io
import sys
from datetime import date, datetime
import requests
import json
import csv
from openpyxl import Workbook
from typing import Any, Tuple
from configparser import ConfigParser
from pathlib import Path


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
        return json.loads(response.content)
    else:
        sys.exit(f"Could not fetch updated records from {url}")


def filter_updated_records(old_file, g2p_data, updated_records_data) -> dict:
    """
    Filter the updated records using the old file submitted to GenCC.
    The relevant updates are those where disease IDs and confidence values have changed.
    Returns a list of records that need to be re-submitted.

    Args:
        old_file (Path): Previous file submitted to GenCC
        g2p_data (list): Data from the G2P all records file
        updated_records_data (dict): Dictionary containing records that have been udpated in G2P since last submission

    Returns:
        dict: Dictionary containing records that need to be re-submitted
    """
    # Prepare output struture
    updated_records_data_filtered = {
        "count": 0,
        "ids": {},
    }

    # Format G2P data for easier access
    g2p_ids_dict = {}
    for row in g2p_data:
        g2p_ids_dict[row["g2p id"]] = row

    # Read file from previous submission and prepare a dictionary for easier access
    old_records = {}
    with open(old_file, mode="r") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            old_records[row["submission_id"]] = row

    for g2p_id in updated_records_data["ids"]:
        submission_id = updated_records_data["ids"][g2p_id]
        if submission_id in old_records:
            old_record = old_records[submission_id]
            # Check disease ID and confidence
            old_disease_id = old_record["disease_id"]
            new_disease_id_mim = g2p_ids_dict[g2p_id]["disease mim"]
            new_disease_id_mondo = g2p_ids_dict[g2p_id]["disease MONDO"]
            old_confidence = old_record["classification_name"]
            new_confidence = g2p_ids_dict[g2p_id]["confidence"]

            if (
                (new_disease_id_mim or new_disease_id_mondo)
                and (
                    old_disease_id != new_disease_id_mim
                    and old_disease_id != new_disease_id_mondo
                )
            ) or old_confidence != new_confidence:
                updated_records_data_filtered["ids"][g2p_id] = submission_id
                updated_records_data_filtered["count"] += 1

    return updated_records_data_filtered


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
        return json.loads(response.content)
    else:
        sys.exit(f"Could not fetch deleted records from {url}")


def post_gencc_submission(list_of_data: list, db_config: dict[str, Any]) -> None:
    """
    Calls the G2P API to write to the db which records are going to be submitted to GenCC.

    Args:
        list_of_data (list): list of records to write to the db
        db_config (dict[str, Any]): dictionary with the database connection details and the API URL

    """
    create_url = f"{db_config['api_url'].rstrip('/')}/gencc_create/"
    login_url = f"{db_config['api_url'].rstrip('/')}/login/"

    login_info = {
        "username": db_config["username"],
        "password": db_config["api_password"],
    }

    response = requests.post(login_url, json=login_info)
    if response.status_code == 200:
        try:
            response_create = requests.post(
                create_url,
                json=list_of_data,
                cookies=response.cookies,
            )
            if response_create.status_code in (200, 201):
                print("GenCC submission created successfully")
            else:
                sys.exit(
                    f"Issues while creating the submissions {response_create.status_code}: {response_create.json}"
                )
        except Exception as e:
            sys.exit("Error:", e)


def write_to_the_GenCC_file(
    records_data: list,
    outfile: Path,
    output_file_issues: Path,
) -> Tuple[Path, list]:
    """
    Write the records to be submitted to the output file 'G2P_GenCC.txt'.
    Records without disease id are excluded from the submission and written
    to file 'record_with_issues.txt'.

    Args:
        records_data (dict[str, Any]): Record of unsubmitted ids
        outfile (Path): Output file to write the submission into
        output_file_issues (Path): File to write the records with issues into

    Returns:
        outfile (Path): The output file with records to be submitted
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
        assertion_criteria_url = (
            "https://www.ebi.ac.uk/gene2phenotype/about/terminology"
        )

        for record in records_data:
            g2p_id = record["g2p id"]

            # Prepare submission ID
            # Create a new submission ID if it's a new record
            if (
                "type_of_submission" not in record
                or record["type_of_submission"] == "create"
            ):
                submission_id = submission_id_base + str(g2p_id[3:])
                type_of_submission = "create"  # new records don't always have the key 'type_of_submission'
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
                issues_with_record[g2p_id] = (
                    f"Unsupported allelic requirement '{record['allelic requirement']}'"
                )
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
        with open(output_file_issues, mode="a") as textfile:
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


def convert_txt_to_excel(input_file: Path, output_file: Path):
    """Converts txt to Excel file

    Args:
        input_file (Path): Text file to be converted
        output_file (Path): Excel file that will be created
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


def prepare_output_dir(path: str) -> Path:
    """
    Create the output directory where files are going to be written.

    Args:
        path (str): path where the G2P and GenCC files are going to be saved

    Returns:
        Path: output directory
    """
    timestamp = create_datetime_now()
    base_dir = Path(path) if path else Path.cwd()
    gencc_dir = base_dir / timestamp
    gencc_dir.mkdir(parents=True, exist_ok=True)

    return gencc_dir


def handle_existing_submission(
    g2p_data: list,
    output_file: Path,
    output_file_updated: Path,
    output_file_issues: Path,
    output_file_deleted: Path,
    old_file: Path,
    db_config: dict[str, Any],
) -> Tuple[Path, Path, list]:
    """
    Method to generate the files for a new submission to GenCC.
    Some records are new, others might have been updated since last submission
    and others might have been deleted in G2P.

    Args:
        g2p_data (list): Data from the G2P all records file
        output_file (Path): Output text file where the submission will be written into
        output_file_updated (Path): Output text file where the updated records will be written into
        output_file_issues (Path): Output text file where the records with issues will be written into
        output_file_deleted (Path): Output text file where the records that have been deleted in G2P will be written into
        old_file (Path): Previous file submitted to GenCC
        db_config (dict[str, Any]): DB/API configuration
    """
    # Retrieves the records that have been updated since last submission
    # We only need to re-submit these if the disease name/IDs have changed
    updated_records_data = get_updated_records(db_config)
    # Review the updated records using the file from the previous submission ('old_file')
    # We only re-submit records where disease IDs or confidence have changed
    updated_records_data_final = filter_updated_records(
        old_file, g2p_data, updated_records_data
    )

    # Retrieves the records that have been deleted in G2P but were submitted to GenCC
    deleted_records_data = get_deleted_records(db_config)

    # List of records that haven't been submitted yet
    # These records are part of public panels and are live in G2P (not deleted)
    unsubmitted = get_unsubmitted_records(db_config)
    # Retrieves the unsubmitted records data from the G2P downloaded file
    records_to_submit = []
    updated_records_to_submit = []
    gencc_list = []

    for row in g2p_data:
        # Prepare unsubmitted records
        if row["g2p id"] in unsubmitted:
            row["type_of_submission"] = "create"
            row["submission_id"] = None
            records_to_submit.append(row)
        if row["g2p id"] in list(updated_records_data_final["ids"].keys()):
            row["type_of_submission"] = "update"
            row["submission_id"] = updated_records_data_final["ids"][row["g2p id"]]
            updated_records_to_submit.append(row)

    # Write the new G2P records to submit to the GenCC file
    outfile, gencc_list_new = write_to_the_GenCC_file(
        records_to_submit, output_file, output_file_issues
    )

    # Write the updated records to a separate file
    outfile_updated_records, gencc_list_updated = write_to_the_GenCC_file(
        updated_records_to_submit, output_file_updated, output_file_issues
    )

    gencc_list = gencc_list_new + gencc_list_updated

    # Write deleted records to a separate file
    with open(output_file_deleted, mode="w") as rw:
        for deleted_record in deleted_records_data["ids"]:
            rw.write(
                "https://www.ebi.ac.uk/gene2phenotype/lgd/" + deleted_record + "\n"
            )

    return outfile, outfile_updated_records, gencc_list


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
        help="Option to determine it is a new submission to GenCC",
    )
    ap.add_argument(
        "--old_file",
        required=False,
        type=Path,
        help="Previous file submitted to GenCC",
    )
    args = ap.parse_args()

    if args.new_submission and args.old_file:
        ap.error("Cannot use --new_submission and --old_file at the same time.")

    if args.old_file and not os.path.isfile(args.old_file):
        sys.exit(f"Invalid file '{args.old_file}'")

    db_config = read_from_config_file(args.config_file)
    skip_write_to_db = args.skip_write_to_db

    # Prepare output directory
    gencc_dir = prepare_output_dir(args.path)
    # File to save records that cannot be submitted due to issues
    output_file_issues = gencc_dir / "records_with_issues.txt"
    # File to save records that are going to be submitted
    output_file = gencc_dir / "G2P_GenCC.txt"
    # File to save records that are going to be submitted as updates
    output_file_updated = gencc_dir / "G2P_GenCC_updated_records.txt"
    # File to save records that have been deleted in G2P
    output_file_deleted = gencc_dir / "G2P_GenCC_deleted_records.txt"
    # File to save records that are going to be submitted (expected format for GenCC)
    final_output_file = gencc_dir / "G2P_GenCC.xlsx"
    # File to save records that are going to be re-submitted (expected format for GenCC)
    final_output_file_updated = gencc_dir / "G2P_GenCC_updated_records.xlsx"

    print("Getting G2P records from the API")
    g2p_data = fetch_g2p_records(db_config)

    if args.new_submission:
        print("Handling new submission")
        # A new submission means that all records are going to be submitted
        # for the first time
        outfile, gencc_list = write_to_the_GenCC_file(
            g2p_data, output_file, output_file_issues
        )
    else:
        print("Handling existing submission")
        outfile, outfile_updated_records, gencc_list = handle_existing_submission(
            g2p_data,
            output_file,
            output_file_updated,
            output_file_issues,
            output_file_deleted,
            args.old_file,
            db_config,
        )

    print("Converting text file to Excel file")
    convert_txt_to_excel(outfile, final_output_file)
    convert_txt_to_excel(outfile_updated_records, final_output_file_updated)

    if not skip_write_to_db:
        print("Writing submission details to the database")
        post_gencc_submission(gencc_list, db_config)


if __name__ == "__main__":
    main()
