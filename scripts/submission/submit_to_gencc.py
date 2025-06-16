#!/usr/bin/env python3

import subprocess 
import argparse
import os
import sys
# import mysql.connector
# from mysql.connector import Error
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


# def get_ols(disease_name):
#     # only using mondo because mondo gives a close enough match with how we name diseases in G2P
#     statement = "No disease mim"
#     endpoint = 'http://www.ebi.ac.uk/ols/api/search?q='
#     ontology = '&ontology=mondo'
#     url = endpoint + disease_name + ontology
#     result = requests.get(url)
#     if result.status_code == 200:
#         final_result = json.loads(result.content)
#         response = final_result["response"]
#         documents = response["docs"]

#         mondo_id = [i["obo_id"] for i in documents if "obo_id" in i ]
#     else:
#         print("Connecting to the " + endpoint + " failed")
#         sys.exit()
#         pass
        
#     if len(mondo_id) > 0:
#         return mondo_id[0]
#     else:
#         return statement

# def fetch_g2p_attribs(host, port, db, user, password):
#     """
#         Fetchs allelic requirement and mutation consequence attribs from the db.
#     """
#     ar_attribs = {}
#     mc_attribs = {}

#     sql_query_attribs = """ SELECT at.code, a.attrib_id, a.value
#                             FROM attrib a 
#                             LEFT JOIN attrib_type at ON at.attrib_type_id = a.attrib_type_id
#                             WHERE at.code = 'allelic_requirement' or at.code = 'mutation_consequence'
#                         """

#     connection = mysql.connector.connect(host=host,
#                                          database=db,
#                                          user=user,
#                                          port=port,
#                                          password=password)

#     try:
#         if connection.is_connected():
#             cursor = connection.cursor()
#             cursor.execute(sql_query_attribs)
#             attrib_data = cursor.fetchall()
#             for row in attrib_data:
#                 if(row[0] == "allelic_requirement"):
#                     ar_attribs[row[1]] = row[2] # key: attrib id; value: attrib name
#                 else:
#                     mc_attribs[row[1]] = row[2] # key: attrib id; value: attrib name
    
#     except Error as e:
#         print("Error while connecting to MySQL", e)
#     finally:
#         if cursor:
#             cursor.close()
#         if connection and connection.is_connected():
#             connection.close()
    
#     return ar_attribs, mc_attribs

def fetch_g2p_records(data: dict[str, Any]):
    fetch_g2p_records = "panel/all/download/"

    url = data["api_url"] + fetch_g2p_records
    response = requests.get(url)

    if response.status_code == 200:
        return response

    else: 
        sys.exit("Issues downloading the file")

    # output_filename = "check.csv"
    # if response.status_code == 200:
    #     with open(output_filename, "wb") as f:
    #         f.write(response.content)
    #     print(f"Downloading G2P files... done. Saved as '{output_filename}'\n")
    #     return output_filename
    # else: 
    #     sys.exit(f"Issues downloading the file: Status Code {response.status_code}")


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

def write_to_the_GenCC_file(data: dict[str, Any], outfile):
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

    return data


def main():
    ap = argparse.ArgumentParser()
    # ap.add_argument("-p", "--path",
    #                 default='/nfs/production/flicek/ensembl/variation/G2P/GenCC_create/',
    #                 help="Path where the G2P and GenCC files are going to be saved")
    ap.add_argument("--config_file", required=True, help="G2P Configuration file")
    args = ap.parse_args()

    # ensembl_dir = os.environ.get('ENSEMBL_ROOT_DIR')
    db_config = read_from_config_file(args.config_file)

    # if not ensembl_dir or not os.path.exists(f"{ensembl_dir}/ensembl-gene2phenotype/scripts/download_file.sh"):
    #     raise FileNotFoundError("ENSEMBL_ROOT_DIR is not set correctly or the script does not exist")

    print("\nDownloading G2P files...")
    file_data = fetch_g2p_records(db_config)
    read_data = reading_data(file_data)
    unsubmitted = get_unsubmitted_record(db_config)
    common = retrieve_unsubmitted_records(read_data, unsubmitted)
    output_file = "G2P_GenCC.txt"
    write_to_the_GenCC_file(common, output_file)

    # outfile = args.path + "/G2P_GenCC.txt"
    # final_output_file = args.path + "/G2P_GenCC.xlsx"

    # final_data_to_submit = {}
    # submitter_id = "GENCC:000112" # G2P submitter id
    # submitter_name = "TGMI G2P" # G2P submitter name
    # assertion_criteria_url = "https://www.ebi.ac.uk/gene2phenotype/terminology"
    # g2p_url = "https://www.ebi.ac.uk/gene2phenotype/"
    # submission_id_base = 1000112 # ID maintained by us. Any alphanumeric string will be accepted up to 64 characters

    # print(f"Fetching all G2P records...")
    # # Fetch the attribs from G2P db
    # ar_attribs, mc_attribs = fetch_g2p_attribs(host, port, db, user, password)
    # # Allelic requeriment and mutation consequence are type SET
    # # Doing a join is complicated when we have multiple values in the SET
    # # That is why we pre-fetched the attribs with fetch_g2p_attribs()
    # all_g2p_records = fetch_g2p_records(host, port, db, user, password, ar_attribs, mc_attribs)
    # print(f"Fetching all G2P records... done\n")

    # # Output format
    # # submission_id: automatic id
    # # hgnc_id: HGNC gene ID
    # # hgnc_symbol: gene symbol (optional)
    # # disease_id: omim, mondo or orpha id
    # # disease_name: our disease name (optional)
    # # moi_id: allelic requeriment GenCC ID
    # # moi_name: allelic requeriment (optional)
    # # submitter_id: G2P GenCC ID
    # # submitter_name: G2P (optional)
    # # classification_id: confidence GenCC ID
    # # classification_name: confidence (optional)
    # # date: format YYYY/MM/DD
    # # public_report_url: G2P URL for the record (optional)
    # # pmids: Listing of PMIDs that have been used in the submission seperated by comma
    # # assertion_criteria_url: G2P terminology url

    # with open(outfile, mode='w') as output_file:
    #     # Write output header
    #     output_file.write("submission_id\thgnc_id\thgnc_symbol\tdisease_id\tdisease_name\tmoi_id\tmoi_name\tsubmitter_id\tsubmitter_name\tclassification_id\tclassification_name\tdate\tpublic_report_url\tpmids\tassertion_criteria_url\n")

    #     for file in files:
    #         file_path = args.path + "/" + file

    #         with gzip.open(file_path, mode='rt') as gz_file:
    #             csv_reader = csv.DictReader(gz_file)
                
    #             for row in csv_reader:
    #                 gene_symbol = row["gene symbol"]
    #                 disease_name = row["disease name"]
    #                 ar = row["allelic requirement"] # we are going to skip records with multiple allelic requirement
    #                 mc_list = row["mutation consequence"].split(";")
    #                 mc_list.sort()
    #                 mc = ";".join(mc_list)
    #                 key = f"{gene_symbol}---{disease_name}---{ar}---{mc}"

    #                 if key not in final_data_to_submit:
    #                     # HGNC ID
    #                     h = row["hgnc id"]
    #                     hgnc_id = f"HGNC:{h}"

    #                     # Prepare MOI - allelic requeriment
    #                     if ar in allelic_requirement:
    #                         moi_id = allelic_requirement[ar]
    #                     else:
    #                         print(f"SKIP: invalid allelic requeriment '{ar}' not found for key '{key}' in file '{file}'")
    #                         continue

    #                     # Get the G2P internal ID to use in the URL
    #                     if key in all_g2p_records:
    #                         g2p_id = all_g2p_records[key]
    #                     else:
    #                         sys.exit(f"Key: '{key}' from download files not found in G2P database")

    #                     record_url = g2p_url + str(g2p_id)

    #                     # the submission id is generated and maintained by us
    #                     g2p_id_formatted = f"{int(g2p_id):05}"
    #                     submission_id = f"{submission_id_base}{g2p_id_formatted}"

    #                     # Confidence
    #                     confidence = row["confidence category"]
    #                     classification_id = confidence_category[confidence]

    #                     # Update disease name to dyadic
    #                     if not disease_name.startswith(gene_symbol):
    #                         new_disease_name = gene_symbol + "-related " + disease_name
    #                     else:
    #                         new_disease_name = disease_name

    #                     # We use the OMIM, Mondo or Orphanet as the disease_id
    #                     disease_id = row["disease mim"]
    #                     disease_ontology = row["disease ontology"]
    #                     if disease_id == "No disease mim" and disease_ontology:
    #                         disease_id = disease_ontology
                        
    #                     # pmids = row["pmids"].replace(";", ",")
    #                     pmids = row["pmids"]

    #                     g2p_date = row["gene disease pair entry date"]
    #                     if g2p_date:
    #                         g2p_date_tmp = datetime.strptime(g2p_date, "%Y-%m-%d %H:%M:%S")
    #                         g2p_date = g2p_date_tmp.strftime("%Y/%m/%d")

    #                     line_to_output = f"{submission_id}\t{hgnc_id}\t{gene_symbol}\t{disease_id}\t{new_disease_name}\t{moi_id}\t{ar}\t{submitter_id}\t{submitter_name}\t{classification_id}\t{confidence}\t{g2p_date}\t{record_url}\t{pmids}\t{assertion_criteria_url}\n"
    #                     output_file.write(line_to_output)

    #                     final_data_to_submit[key] = 1

    # output_file.close()

    # convert_txt_to_excel(outfile, final_output_file)

if __name__ == '__main__':
    main()