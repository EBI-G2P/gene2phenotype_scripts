#!/usr/bin/env python3

import subprocess 
import argparse
import os
import sys
from datetime import date, datetime
import requests
# import response
import json
import gzip
import csv
# from openpyxl import Workbook
from typing import Any
from configparser import ConfigParser
import io
# Mapping terms to GenCC IDs
allelic_requirement = {
    "biallelic_autosomal" : "HP:0000007",
    "monoallelic_autosomal" : "HP:0000006",
    "monoallelic_X_hemizygous"  : "HP:0001417",
    "monoallelic_Y_hemizygous" : "HP:0001450",
    "monoallelic_X_heterozygous" : "HP:0001417",
    "mitochondrial" : "HP:0001427",
    "monoallelic_PAR" : "HP:0000006",
    "biallelic_PAR" : "HP:0000007"
}
confidence_category = {
    "definitive" : "GENCC:100001",
    "strong" : "GENCC:100002",
    "moderate" : "GENCC:100003",
    "limited" : "GENCC:100004",
    "disputed": "GENCC:100005",
    "refuted": "GENCC:100006"
}


def fetch_g2p_records(data: dict[str, Any]):
    fetch_g2p_records = "panel/all/download/"

    url = data["api_url"] + fetch_g2p_records
    response = requests.get(url)

    if response.status_code == 200:
        return response

    else: 
        sys.exit("Issues downloading the file")

def reading_data(response) -> list:
    csv_content = io.StringIO(response.content.decode("utf-8"))
    reader = csv.DictReader(csv_content)

    return reader

def get_unsubmitted_record(data: dict[str, Any])-> list:
    fetch_unsubmitted_record = "unsubmitted-stable-ids/"

    url = data["api_url"] + fetch_unsubmitted_record
    response = requests.get(url)
    if response.status_code == 200:
        unsubmitted = json.loads(response.content)
        unsubmitted_ids = {record["stable_id"] for record in unsubmitted}
        return unsubmitted_ids
    else:
        print("Could not fetch unsubmitted records")

def retrieve_unsubmitted_records(data: dict[str, Any], unsubmitted: list) -> list:
    records = [row for row in data if row["g2p id"] in unsubmitted]

    return records

def create_gencc_submission_record(data:dict[str, Any], submission_id, date, type_of, stable_id):
    create_url = "gencc_create/"

    login_url = "login/"

    login_info = {"username": data["username"], "password": data["password"]}

    response = requests.post(data["api_url"]+login_url, json=login_info)
    create_info = {
        "submission_id": submission_id,
        "date_of_submission": date,
        "type_of_submission": type_of,
        "g2p_stable_id" : stable_id
    }
    if response.status_code == 200:
        try:
            response_create = requests.post(
                data["api_url"]+create_url, json=create_info, cookies=response.cookies
            )
            if response_create.status_code in (200, 201):
                print(f"GenCC submission for the record {stable_id} was created successfully")
            else:
                print(f"Issues creating the submissions {response_create.status_code} {response_create.json}")
        except Exception as e:
            print('Error:', e)

def write_to_the_GenCC_file(data: dict[str, Any], outfile, dry, db_config):
    with open(outfile, mode='w') as output_file:
         output_file.write("submission_id\thgnc_id\thgnc_symbol\tdisease_id\tdisease_name\tmoi_id\tmoi_name\tsubmitter_id\tsubmitter_name\tclassification_id\tclassification_name\tdate\tpublic_report_url\tpmids\tassertion_criteria_url\n")
        
         submission_id_base = "1000112"
         for record in data:
             g2p_id = record['g2p id']
             submission_id = submission_id_base + str(g2p_id[3:])

             hgnc_id = record['hgnc id']
             hgnc_symbol = record['gene symbol']

             disease_id = record['disease mim'] or record['disease MONDO']
             disease_name = record['disease name']
             moi_id = allelic_requirement[record['allelic requirement']]
             moi_name = record['allelic requirement']
             submitter_id = 'GENCC:000112'
             submitter_name = 'TGMI G2P'
             classification_id = confidence_category[record['confidence']]
             classification_name = record["confidence"]
             dt = datetime.fromisoformat(record['date of last review'])
             date = dt.strftime("%Y/%m/%d")
             record_url = "https://www.ebi.ac.uk/gene2phenotype/lgd/" + g2p_id
             pmids = record['publications']
             assertion_criteria_url = "https://www.ebi.ac.uk/gene2phenotype/terminology"

             line_to_output = f"{submission_id}\t{hgnc_id}\t{hgnc_symbol}\t{disease_id}\t{disease_name}\t{moi_id}\t{moi_name}\t{submitter_id}\t{submitter_name}\t{classification_id}\t{classification_name}\t{date}\t{record_url}\t{pmids}\t{assertion_criteria_url}\n"
             output_file.write(line_to_output)
             if dry == "False":
                 type_of = "create"
                 db_date = dt.strftime("%Y-%m-%d") 
                 create_gencc_submission_record(db_config, submission_id, db_date, type_of, g2p_id)
    return output_file



def convert_txt_to_excel(input_file, output_file):
    """
        Converts a text file to an Excel file.
    """
    wb = Workbook()
    ws = wb.active

    delimiter='\t'

    with open(input_file, 'r') as file:
        for row_index, line in enumerate(file, start=1):
            for col_index, value in enumerate(line.strip().split(delimiter), start=1):
                ws.cell(row=row_index, column=col_index, value=value)
    wb.save(output_file)

def read_from_config_file(config_file: str) -> dict[str, Any]:
    """_summary_

    Args:
        config_file (str): _description_

    Returns:
        dict[str, Any]: _description_
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
    data["password"] = config["api"]["password"]

    return data


def main():
    ap = argparse.ArgumentParser()
    # ap.add_argument("-p", "--path",
    #                 default='/nfs/production/flicek/ensembl/variation/G2P/GenCC_create/',
    #                 help="Path where the G2P and GenCC files are going to be saved")
    ap.add_argument("--config_file", required=True, help="G2P Configuration file")
    ap.add_argument("--dry", required=False, help="If dry is False, it creates an actual GenCC submission so updates the db" )
    args = ap.parse_args()

    # ensembl_dir = os.environ.get('ENSEMBL_ROOT_DIR')
    db_config = read_from_config_file(args.config_file)
   


    if args.dry:
        dry = args.dry

    # if not ensembl_dir or not os.path.exists(f"{ensembl_dir}/ensembl-gene2phenotype/scripts/download_file.sh"):
    #     raise FileNotFoundError("ENSEMBL_ROOT_DIR is not set correctly or the script does not exist")

    print("\nDownloading G2P files...")
    file_data = fetch_g2p_records(db_config)
    read_data = reading_data(file_data)
    unsubmitted = get_unsubmitted_record(db_config)
    common = retrieve_unsubmitted_records(read_data, unsubmitted)
    output_file = "G2P_GenCC.txt"
    outfile = write_to_the_GenCC_file(common, output_file, dry, db_config)
    #lonvert_txt_to_excel(outfile, final_output_file)

if __name__ == '__main__':
    main()