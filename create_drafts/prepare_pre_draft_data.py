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


"""
Script to read data from ClinGen or PanelApp and create JSON file with record data
that is going to be used to create the draft records.

For each ClinGen or PanelApp record, the script checks if the gene-disease association
is already in G2P and if not, it adds it to the list of records to create as draft records.
If the gene-disease association is already in G2P but there are differences in the disease
or mechanism evidence, it adds it to a list of potential records to create as draft records
that need to be checked.

Options:
        --source         ClinGen or PanelApp
        --input_file     Data to be used to create drafts (supported formats: json and tsv)
                         The ClinGen file is the output file of the script gemini_analise_clingen.py (json)
                         The PanelApp file can be downloaded from https://panelapp.genomicsengland.co.uk/panels/112/ (tsv)

Example how to run:
    python prepare_pre_draft_data.py --source ClinGen --input_file clingen_extracted_data.json
"""


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


def download_g2p_records(
    api_url: str, cookies: requests.cookies.RequestsCookieJar
) -> list:
    url = f"{api_url.rstrip('/')}/panel/all/download/"

    try:
        response = requests.get(url, cookies=cookies)
    except Exception as e:
        sys.exit("Failed to download G2P file:", str(e))
    else:
        if response.status_code == 200:
            csv_content = io.StringIO(response.content.decode("utf-8"))
            reader = csv.DictReader(csv_content)

            return list(reader)
        else:
            sys.exit(
                f"Failed to download G2P file. Status code: {response.status_code}"
            )


def process_g2p_records(g2p_records) -> list:
    """
    Format the g2p records data
    """
    g2p_data_by_gene = {}

    for record in g2p_records:
        if record["gene symbol"] not in g2p_data_by_gene:
            g2p_data_by_gene[record["gene symbol"]] = [
                {
                    "g2p_id": record["g2p id"],
                    "disease": record["disease name"],
                    "mondo_id": record["disease MONDO"],
                    "omim_id": record["disease mim"],
                    "panel": record["panel"],
                }
            ]
        else:
            g2p_data_by_gene[record["gene symbol"]].append(
                {
                    "g2p_id": record["g2p id"],
                    "disease": record["disease name"],
                    "mondo_id": record["disease MONDO"],
                    "omim_id": record["disease mim"],
                    "panel": record["panel"],
                }
            )

    return g2p_data_by_gene


def process_clingen_data(file: str, g2p_data_by_gene: list, output_file: str) -> None:
    """
    Read the ClinGen json file
    """
    total_records = 0
    new_records = 0
    records_to_check = 0
    drafts_to_create = []

    with open(file) as fh:
        data = json.load(fh)
        for record in data:
            found_disease = False
            if record["clingen_panel"] != "Hearing Loss Gene Curation Expert Panel":
                continue
            total_records += 1

            # Format the mechanism evidence to the expected format
            new_mechanism_evidence = []
            for mechanism_evidence in record["evidence"]:
                if isinstance(mechanism_evidence, dict):
                    if "description" in mechanism_evidence:
                        new_string = f"{mechanism_evidence['type']}: {mechanism_evidence['description']}"
                    elif "evidence" in mechanism_evidence:
                        new_string = f"{mechanism_evidence['type']}: {mechanism_evidence['evidence']}"
                    else:
                        new_string = f"{mechanism_evidence['type']}"
                    new_mechanism_evidence.append(new_string)
            if new_mechanism_evidence:
                record["evidence"] = new_mechanism_evidence

            gene = record["gene_symbol"]
            # The gene is already associated with g2p records
            if gene in g2p_data_by_gene:
                for g2p_record in g2p_data_by_gene[gene]:
                    # Check if it is the same disease using the mondo ids
                    if g2p_record["mondo_id"] == record["mondo_id"]:
                        found_disease = True
                    else:
                        # Check if it is the same disease
                        new_g2p_disease = re.sub(
                            ".*\-related ", "", g2p_record["disease"]
                        )
                        new_clingen_disease = record["disease"].replace(",", "")
                        new_clingen_disease = re.sub(
                            " type \d+", "", new_clingen_disease
                        )
                        if new_g2p_disease.strip() == new_clingen_disease.strip():
                            found_disease = True
                if not found_disease:
                    records_to_check += 1
                    record["type"] = "probably_create"
                    record["g2p_data"] = g2p_data_by_gene[gene]
                    drafts_to_create.append(record)
            else:
                # Gene not yet curated in G2P
                new_records += 1
                record["type"] = "create"
                record["g2p_data"] = "NA"
                drafts_to_create.append(record)

    print(f"\nTotal records: {total_records}")
    print(f"New records: {new_records}")
    print(f"Potential new records: {records_to_check}")

    with open(output_file, "wt") as wr:
        json.dump(drafts_to_create, wr, indent=2)


def normalize_disease_name_for_comparison(disease_name: str) -> str:
    """
    Normalize disease names before comparing G2P and PanelApp records.
    """
    disease_name = re.sub(".*\\-related ", "", disease_name)
    disease_name = disease_name.replace(",", "")
    disease_name = re.sub(r"\s+type\s+\d+[A-Za-z]*\b", "", disease_name)
    disease_name = re.sub(r"\s+\d+[A-Za-z]+\W*$", "", disease_name)
    disease_name = re.sub(r"\s+dominant\W*$", "", disease_name, flags=re.IGNORECASE)
    disease_name = re.sub(r"\s+", " ", disease_name)
    return disease_name.strip().lower()


def disease_word_order_key(disease_name: str) -> str:
    """
    Return a comparison key that ignores word order.
    """
    words = re.findall(r"[a-z0-9]+", disease_name)
    return " ".join(sorted(words))


def panelapp_disease_matches_g2p(g2p_disease: str, panelapp_diseases: list) -> bool:
    """
    Compare a G2P disease name with PanelApp disease names.
    """
    new_g2p_disease = normalize_disease_name_for_comparison(g2p_disease)
    for disease in panelapp_diseases:
        new_panelapp_disease = normalize_disease_name_for_comparison(disease["name"])

        if (
            new_g2p_disease == new_panelapp_disease
            or disease_word_order_key(new_g2p_disease)
            == disease_word_order_key(new_panelapp_disease)
        ):
            return True

    return False


def process_panelapp_data(
    file: str,
    g2p_data_by_gene: list,
    output_file: str,
    clingen_file: str,
    clingen_panel: str,
    katherine_list: dict,
) -> None:
    total_records = 0
    new_records = 0
    drafts_to_create = []
    clingen_records = {}

    # Load the ClinGen data
    with open(clingen_file) as fh:
        data = json.load(fh)
        for record in data:
            if record["gene_symbol"] in clingen_records:
                clingen_records[record["gene_symbol"]].append(record)
            else:
                clingen_records[record["gene_symbol"]] = [record]

    with open(file) as fh:
        data = csv.DictReader(fh, delimiter="\t", quotechar='"')
        for row in data:
            if row["Entity type"] != "gene":
                continue
            total_records += 1
            gene_symbol = row["Gene Symbol"]
            hgnc_id = row["HGNC"]
            omim_id = row["Omim"]  # gene id
            model_of_inheritance = row["Model_Of_Inheritance"]
            phenotypes = row["Phenotypes"]
            publications = row["Publications"]
            disease_ids = []

            # Format phenotypes
            list_phenotypes = phenotypes.split(";")
            phenotypes_final = []
            for pheno in list_phenotypes:
                pheno_search_omim = re.findall(
                    r"(OMIM:\s*\d{6})|(MIM#\d{6})|(\s\d{6}$)", pheno
                )
                pheno = re.sub(r",*\s*OMIM:\s*\d{6}", "", pheno)
                pheno = re.sub(r",*\s*MIM#\d{6}", "", pheno)
                pheno = re.sub(r",*\s\d{6}$", "", pheno)
                for pheno_omim_id in pheno_search_omim:
                    for idx, element in enumerate(pheno_omim_id):
                        if element:
                            omim_id_clean = element.strip().replace("OMIM:", "")
                            phenotypes_final.append(
                                {"name": pheno, "disease_id": omim_id_clean.strip()}
                            )
                            # Add disease id to list
                            disease_ids.append(omim_id_clean.strip())
                # Add the phenotypes that do not have omim id
                if not pheno_search_omim:
                    phenotypes_final.append({"name": pheno})

            # Format publications
            list_publications = publications.split(";")
            publications_final = []
            for publication in list_publications:
                if (
                    publication
                    == "Ostergaard et al., 2015, J. Med. Genet., 52, 203-207."
                ):
                    publication = "25604084"
                elif (
                    publication
                    == "Reyes et al., 2005, Am. J. Hum. Genet., 97,  186-193."
                ):
                    publication = "29572490"

                publication = re.sub(r"PMID:\s*", "", publication)
                publication = re.sub(r"PubMed:\s*", "", publication)
                publication = publication.replace(
                    "Besse et al., 2015, Cell Metab., 21(3), 417-427", ""
                )
                if publication and not publication[0].isalpha():
                    publications_final.append(publication.strip())

            draft_record = {
                "gene_symbol": gene_symbol,
                "hgnc_id": hgnc_id,
                "omim_id": omim_id,
                "allelic_requirement": model_of_inheritance,
                "pmids": publications_final,
                "phenotypes": phenotypes_final,
                "disease_id": disease_ids,
            }

            # Add ClinGen data to the record
            clingen_info = []
            if gene_symbol in clingen_records:
                clingen_list = clingen_records[gene_symbol]
                for clingen_record in clingen_list:
                    if clingen_record["clingen_panel"] == clingen_panel:
                        clingen_info.append(clingen_record)
            if clingen_info:
                draft_record["clingen_data"] = clingen_info
            if katherine_list and gene_symbol in katherine_list:
                draft_record["katherine_schon_data"] = katherine_list[gene_symbol]

            # Check if gene is curated in G2P
            draft_record["g2p_data"] = []
            found = False
            if gene_symbol in g2p_data_by_gene:
                for g2p_record in g2p_data_by_gene[gene_symbol]:
                    # Match records using the disease ids (omim or mondo)
                    if (g2p_record["omim_id"] in draft_record["disease_id"]
                        or g2p_record["mondo_id"] in draft_record["disease_id"]):
                        found = True
                    # Match records using the disease name
                    elif panelapp_disease_matches_g2p(
                        g2p_record["disease"], draft_record["phenotypes"]
                    ):
                        found = True
                    draft_record["g2p_data"].append(g2p_record)
                if not found:
                    draft_record["type"] = "probably_create"
                else:
                    draft_record["type"] = "in_g2p"
            else:
                draft_record["type"] = "create"

            drafts_to_create.append(draft_record)
            new_records += 1

        with open(output_file, "wt") as wr:
            json.dump(drafts_to_create, wr, indent=2)


def read_extra_file(katherine_file: str) -> dict:
    """
    Read Katherine Schon MT gene list
    """
    final_list = {}

    if not katherine_file:
        return

    if not os.path.isfile(katherine_file):
        sys.exit(f"Invalid Katherine Schon MT gene list file '{katherine_file}'")

    with open(katherine_file, "r", encoding="utf-8") as fh:
        for row in fh:
            data_row = row.split("\t")
            if data_row[0] not in final_list:
                final_list[data_row[0]] = [
                    {
                        "pmids": data_row[1].strip().replace('"', ""),
                        "comment": data_row[2].strip().replace('"', ""),
                        "gene_rating": data_row[3].strip().replace('"', ""),
                    }
                ]
            else:
                final_list[data_row[0]].append(
                    {
                        "pmids": data_row[1].strip().replace('"', ""),
                        "comment": data_row[2].strip().replace('"', ""),
                        "gene_rating": data_row[3].strip().replace('"', ""),
                    }
                )

    return final_list


def main():
    parser = argparse.ArgumentParser(description="Create draft records")
    parser.add_argument("--source", required=True, help="ClinGen or PanelApp")
    parser.add_argument(
        "--input_file",
        required=True,
        help="Data to be used to create drafts (supported formats: json and tsv)",
    )
    parser.add_argument(
        "--clingen_file", required=True, help="ClinGen data to append to PanelApp records (supported format: json)"
    )
    parser.add_argument(
        "--clingen_panel", required=True, help="ClinGen panel to append to PanelApp records"
    )
    parser.add_argument(
        "--katherine_file",
        required=False,
        help="Katherine Schon MT gene list (supported format: json)",
    )
    parser.add_argument(
        "--output_file",
        required=False,
        help="Output file name for pre-draft records (supported format: json)",
    )
    parser.add_argument(
        "--config", required=True, help="Config file with details to G2P API"
    )
    args = parser.parse_args()

    input_file = args.input_file
    source = args.source
    clingen_file = args.clingen_file
    clingen_panel = args.clingen_panel
    katherine_file = args.katherine_file
    output_file = args.output_file

    if not output_file:
        output_file = f"{source}_pre_draft_records.json"

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

    # Download latest G2P data
    cookies = login(api_username, api_password, api_url)
    g2p_records = download_g2p_records(api_url, cookies)
    g2p_data_by_gene = process_g2p_records(g2p_records)
    logout(api_url, cookies)

    if source.lower() == "clingen":
        process_clingen_data(input_file, g2p_data_by_gene, output_file)

    if source.lower() == "panelapp":
        katherine_list = {}
        if katherine_file:
            katherine_list = read_extra_file(katherine_file)
        process_panelapp_data(
            input_file, g2p_data_by_gene, output_file, clingen_file, clingen_panel, katherine_list
        )


if __name__ == "__main__":
    main()
