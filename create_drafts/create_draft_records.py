#!/usr/bin/env python3

import sys
import argparse
import configparser
import csv
import io
import os.path
import re
import requests
import json


"""
Script to create automatic draft records in G2P.

For ClinGen, the input file must be the JSON output from
gemini_analise_clingen.py. Records from fetch_clingen.py are not sufficient
because they do not contain the Gemini-extracted publications, mechanism,
allelic requirement, phenotypes, and mechanism evidence required for draft
creation.

Supported sources:
- clingen creates drafts from Gemini-enriched ClinGen JSON.
- immuno creates drafts from an immunology workbook plus Gemini-enriched ClinGen JSON.
- inhouse creates drafts from an in-house mined-publication CSV.

Options:
        --source             clingen, immuno or inhouse
        --input_file         ClinGen JSON or inhouse CSV
        --config             Config file with G2P API details
        --panel              G2P panel name e.g. "Ear disorders"
        --clingen_panel      Optional ClinGen panel filter when --source is clingen
        --immuno_file        Immunology workbook when --source is immuno
        --genes_to_include   Optional gene-symbol allowlist, one gene per line
        --output_file        Optional JSON output file; if set, drafts are not inserted into the API

Example ClinGen usage:
    python create_draft_records.py \
        --config config.ini \
        --source clingen \
        --panel "Ear disorders" \
        --clingen_panel "Hearing Loss Gene Curation Expert Panel" \
        --input_file clingen_extracted_data_2026-09-25_gemini.json

Example Immuno usage:
    python create_draft_records.py \
        --config config.ini \
        --source immuno \
        --panel "Immunological disorders" \
        --input_file clingen_extracted_data_2026-09-25_gemini.json \
        --immuno_file immuno_panel.xlsx

Example inhouse usage:
    python create_draft_records.py \
        --config config.ini \
        --source inhouse \
        --panel "Cardiac disorders" \
        --input_file mined_publications.csv
"""

valid_mechanisms = [
    "loss of function",
    "gain of function",
    "dominant negative",
    "undetermined non loss of function",
    "undetermined",
]

allelic_requirement_mapping = {
    "autosomal recessive": "biallelic_autosomal",
    "biallelic": "biallelic_autosomal",
    "autosomal dominant": "monoallelic_autosomal",
    "X-linked": "monoallelic_X",
}

MAX_SESSION_NAME_LENGTH = 100
REQUIRED_GEMINI_CLINGEN_FIELDS = {
    "pmids",
    "disease_id",
    "mechanism",
    "allelic_requirement",
    "phenotypes",
    "evidence",
    "comment",
}


def build_session_name(*parts: str) -> str:
    """Build a DB-safe session name within the column length limit."""
    return "_".join(parts)[:MAX_SESSION_NAME_LENGTH]


def normalize_excel_headers(headers: tuple) -> list[str]:
    """Make worksheet headers usable as unique dictionary keys."""
    counts = {}
    normalized_headers = []

    for index, header in enumerate(headers, start=1):
        header_text = str(header).strip() if header is not None else f"column_{index}"
        if not header_text:
            header_text = f"column_{index}"

        counts[header_text] = counts.get(header_text, 0) + 1
        if counts[header_text] > 1:
            header_text = f"{header_text}_{counts[header_text]}"

        normalized_headers.append(header_text)

    return normalized_headers


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
    except Exception as error:
        sys.exit(f"Failed to download G2P file: {error}")

    if response.status_code != 200:
        sys.exit(f"Failed to download G2P file. Status code: {response.status_code}")

    csv_content = io.StringIO(response.content.decode("utf-8"))
    reader = csv.DictReader(csv_content)
    return list(reader)


def process_g2p_records(g2p_records: list) -> dict:
    """Index G2P records by gene symbol for ClinGen duplicate checks."""
    g2p_data_by_gene = {}

    for record in g2p_records:
        gene_symbol = record["gene symbol"]
        g2p_data_by_gene.setdefault(gene_symbol, []).append(
            {
                "g2p_id": record["g2p id"],
                "disease": record["disease name"],
                "mondo_id": record["disease MONDO"],
                "omim_id": record["disease mim"],
                "panel": record["panel"],
            }
        )

    return g2p_data_by_gene


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


def handle_draft_record(
    final_draft: dict,
    api_url: str,
    cookies: requests.cookies.RequestsCookieJar,
    output_records: list | None,
) -> None:
    draft_data = {"json_data": final_draft, "status": "automatic"}
    if output_records is not None:
        output_records.append(draft_data)
    else:
        insert_draft_record(api_url, cookies, draft_data)


def write_draft_output(output_file: str, output_records: list) -> None:
    with open(output_file, "wt", encoding="utf-8") as fh:
        json.dump(output_records, fh, indent=2)


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


def normalize_disease_name_for_comparison(disease_name: str) -> str:
    """Normalize disease names before comparing G2P and ClinGen records."""
    disease_name = re.sub(".*\\-related ", "", disease_name)
    disease_name = disease_name.replace(",", "")
    disease_name = re.sub(r"\s+type\s+\d+[A-Za-z]*\b", "", disease_name)
    disease_name = re.sub(r"\s+\d+[A-Za-z]+\W*$", "", disease_name)
    disease_name = re.sub(r"\s+dominant\W*$", "", disease_name, flags=re.IGNORECASE)
    disease_name = re.sub(r"\s+", " ", disease_name)
    return disease_name.strip().lower()


def disease_word_order_key(disease_name: str) -> str:
    """Return a comparison key that ignores word order."""
    words = re.findall(r"[a-z0-9]+", disease_name)
    return " ".join(sorted(words))


def disease_matches_g2p(g2p_disease: str, source_disease: str) -> bool:
    normalized_g2p_disease = normalize_disease_name_for_comparison(g2p_disease)
    normalized_source_disease = normalize_disease_name_for_comparison(source_disease)

    return (
        normalized_g2p_disease == normalized_source_disease
        or disease_word_order_key(normalized_g2p_disease)
        == disease_word_order_key(normalized_source_disease)
    )


def read_genes_to_include(genes_to_include: str) -> list:
    """Read file with list of genes to include in draft records."""
    if not os.path.isfile(genes_to_include):
        sys.exit(f"Invalid file with genes to include '{genes_to_include}'")

    with open(genes_to_include, "r", encoding="utf-8") as fh:
        return [line.strip() for line in fh if line.strip()]


def format_mechanism_evidence(record: dict) -> list:
    formatted_evidence = []
    for mechanism_evidence in record.get("evidence", []):
        if isinstance(mechanism_evidence, dict):
            mechanism_evidence_type = (
                mechanism_evidence.get("type", "")
                or mechanism_evidence.get("evidence_type", "")
            )
            if "description" in mechanism_evidence:
                formatted_evidence.append(
                    f"{mechanism_evidence_type}: {mechanism_evidence['description']}"
                )
            elif "evidence" in mechanism_evidence:
                formatted_evidence.append(
                    f"{mechanism_evidence_type}: {mechanism_evidence['evidence']}"
                )
            else:
                formatted_evidence.append(mechanism_evidence_type)
        else:
            formatted_evidence.append(str(mechanism_evidence))

    return formatted_evidence


def build_clingen_mechanism_evidence(record: dict) -> str:
    evidence_parts = []
    mechanism_evidence = format_mechanism_evidence(record)

    if mechanism_evidence:
        evidence_parts.append("Summary: " + "; ".join(mechanism_evidence))

    for exp_evidence in record.get("experimental_evidence", []):
        pmid_match = re.search(r"PMID:\s*(\d+)", exp_evidence.get("reference", ""))
        pmid = pmid_match.group(1) if pmid_match else None
        if not pmid:
            continue

        evidence_parts.append(
            f"Experimental evidence for PMID {pmid}: "
            f"{exp_evidence.get('experimental_category', '')}; "
            f"Explanation: {exp_evidence.get('explanation', '')}"
        )

    return "\n\n".join(evidence_parts)


def validate_gemini_clingen_record(record: dict) -> None:
    missing_fields = REQUIRED_GEMINI_CLINGEN_FIELDS - set(record)
    if missing_fields:
        missing_fields_list = ", ".join(sorted(missing_fields))
        gene_symbol = record.get("gene_symbol", "unknown gene")
        sys.exit(
            "ERROR: ClinGen input must be generated by gemini_analise_clingen.py. "
            f"Record for {gene_symbol} is missing fields: {missing_fields_list}"
        )


def validate_clingen_panel(records: list, clingen_panel: str) -> None:
    if not clingen_panel:
        return

    available_panels = sorted(
        {record.get("clingen_panel", "") for record in records if record.get("clingen_panel")}
    )
    if clingen_panel not in available_panels:
        available_panel_text = "\n".join(f"- {panel}" for panel in available_panels)
        sys.exit(
            f"ERROR: ClinGen panel '{clingen_panel}' was not found in input.\n"
            f"Available panels:\n{available_panel_text}"
        )


def matching_g2p_records(record: dict, g2p_data_by_gene: dict) -> list:
    matched_records = []
    gene = record["gene_symbol"]
    record_mondo_id = record.get("mondo_id", "")

    for g2p_record in g2p_data_by_gene.get(gene, []):
        g2p_mondo_id = g2p_record.get("mondo_id", "")
        if record_mondo_id and g2p_mondo_id and g2p_mondo_id == record_mondo_id:
            matched_records.append(g2p_record)
        elif disease_matches_g2p(g2p_record["disease"], record["disease"]):
            matched_records.append(g2p_record)

    return matched_records


def build_existing_g2p_comment(g2p_records: list) -> str:
    existing_records = []
    for record in g2p_records:
        existing_records.append(
            f"{record['g2p_id']} | {record['panel']} | {record['disease']}"
        )
    return "Existing G2P records for this gene: " + "; ".join(existing_records) + "\n"


def build_disease_cross_references(record: dict, api_url: str) -> list:
    disease_cross_references = []

    if record.get("mondo_id"):
        disease_data = fetch_disease_info(api_url, record["mondo_id"])
        disease_cross_references.append(
            {
                "source": "Mondo",
                "identifier": record["mondo_id"],
                "disease_name": disease_data["disease"].lower(),
                "original_disease_name": disease_data["disease"],
            }
        )

    disease_ids = record.get("disease_id", [])
    if isinstance(disease_ids, str):
        disease_ids = disease_ids.split(",")

    for disease_id in disease_ids:
        disease_id = str(disease_id).strip()
        if not disease_id:
            continue

        if disease_id.startswith("OMIM:") or disease_id.startswith("MIM:"):
            disease_id = disease_id.replace("OMIM:", "").replace("MIM:", "").strip()
            if not disease_id.isdigit():
                continue
            disease_data = fetch_disease_info(api_url, disease_id)
            disease_cross_references.append(
                {
                    "source": "OMIM",
                    "identifier": disease_id,
                    "disease_name": disease_data["disease"].lower(),
                    "original_disease_name": disease_data["disease"],
                }
            )
        elif disease_id.startswith("MONDO:"):
            disease_data = fetch_disease_info(api_url, disease_id)
            disease_cross_references.append(
                {
                    "source": "Mondo",
                    "identifier": disease_id,
                    "disease_name": disease_data["disease"].lower(),
                    "original_disease_name": disease_data["disease"],
                }
            )

    return disease_cross_references


def prepare_clingen_draft_records(
    input_file: str,
    panel_name: str,
    source: str,
    automatic_drafts: dict,
    api_url: str,
    cookies: requests.cookies.RequestsCookieJar,
    clingen_panel: str,
    genes_to_include: str,
    output_records: list | None = None,
) -> None:
    with open(input_file) as fh:
        data = json.load(fh)

    validate_clingen_panel(data, clingen_panel)
    for record in data:
        if not clingen_panel or record.get("clingen_panel") == clingen_panel:
            validate_gemini_clingen_record(record)

    g2p_records = download_g2p_records(api_url, cookies)
    g2p_data_by_gene = process_g2p_records(g2p_records)
    genes_to_keep = read_genes_to_include(genes_to_include) if genes_to_include else []
    genes_found = []
    total_records = 0
    new_records = 0
    records_to_check = 0
    skipped_existing_records = 0

    for record in data:
        if clingen_panel and record.get("clingen_panel") != clingen_panel:
            continue

        total_records += 1
        gene = record["gene_symbol"]

        if genes_to_keep and gene not in genes_to_keep:
            continue

        genes_found.append(gene)
        matched_g2p_records = matching_g2p_records(record, g2p_data_by_gene)
        if matched_g2p_records:
            skipped_existing_records += 1
            continue

        existing_gene_records = g2p_data_by_gene.get(gene, [])
        if existing_gene_records:
            records_to_check += 1
        else:
            new_records += 1

        final_draft = {}
        final_draft["extra_comment"] = ""
        if existing_gene_records:
            final_draft["extra_comment"] += build_existing_g2p_comment(
                existing_gene_records
            )
        if record.get("comment"):
            final_draft["extra_comment"] += f"Gemini comment: {record['comment']}\n"

        final_draft["locus"] = gene
        final_draft["panels"] = [panel_name]
        final_draft["public_comment"] = ""
        final_draft["private_comment"] = ""
        final_draft["cross_cutting_modifier"] = []
        final_draft["variant_types"] = []
        final_draft["variant_descriptions"] = []
        final_draft["variant_consequences"] = []

        final_draft["source_data"] = {"name": source}
        if "url" in record:
            final_draft["source_data"]["url"] = record["url"]

        final_draft["confidence"] = ""
        if record.get("confidence"):
            final_draft["source_data"]["confidence"] = record["confidence"].lower()

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
            final_draft["source_data"]["phenotypes"] = record["phenotypes"]

        final_draft["allelic_requirement"] = ""
        if record["allelic_requirement"] in allelic_requirement_mapping:
            final_draft["allelic_requirement"] = allelic_requirement_mapping[
                record["allelic_requirement"]
            ]
        else:
            final_draft["extra_comment"] += (
                "Unsupported allelic requirement: "
                + record["allelic_requirement"]
                + "\n"
            )

        final_draft["source_data"]["mechanism_comment"] = ""
        valid_mechanism = ""
        if (
            record.get("mechanism")
            and record["mechanism"].lower().replace("-", " ") not in valid_mechanisms
        ):
            final_draft["source_data"]["mechanism_comment"] += (
                "Unsupported mechanism: " + record["mechanism"]
            )
        elif record.get("mechanism"):
            valid_mechanism = record["mechanism"].lower().replace("-", " ")

        final_draft["molecular_mechanism"] = {
            "name": "",
            "support": "",
        }
        final_draft["mechanism_synopsis"] = []
        final_draft["mechanism_evidence"] = []
        final_draft["source_data"]["mechanism"] = valid_mechanism

        mechanism_evidence = build_clingen_mechanism_evidence(record)
        if mechanism_evidence:
            final_draft["source_data"]["mechanism_evidence"] = mechanism_evidence

        final_draft["disease"] = {
            "disease_name": "",
            "cross_references": [],
        }
        final_draft["source_data"]["disease"] = record["disease"]
        final_draft["source_data"]["disease_cross_references"] = (
            build_disease_cross_references(record, api_url)
        )

        final_draft["session_name"] = build_session_name(
            gene,
            final_draft["allelic_requirement"],
            record["disease"],
        )

        if final_draft["session_name"] in automatic_drafts:
            print(
                f"Draft record for session '{final_draft['session_name']}' already exists. Skipping insertion."
            )
            continue

        handle_draft_record(
            final_draft,
            api_url,
            cookies,
            output_records,
        )

    print(f"\nTotal ClinGen records in panel: {total_records}")
    print(f"New records: {new_records}")
    print(f"Potential new records: {records_to_check}")
    print(f"Existing records skipped: {skipped_existing_records}")

    if genes_to_keep:
        genes_not_found = set(genes_to_keep) - set(genes_found)
        if genes_not_found:
            print(
                "\nThe following genes were not found in the ClinGen data but are "
                f"in the list of genes to keep:\n{', '.join(sorted(genes_not_found))}"
            )


def prepare_draft_records_pipeline(
    input_file: str,
    g2p_to_pmids: dict,
    panel_name: str,
    source: str,
    automatic_drafts: dict,
    api_url: str,
    cookies: requests.cookies.RequestsCookieJar,
    output_records: list | None = None,
) -> None:
    """
    Read the input csv file. This is the fake G2P file that was used by the pipeline.
    """

    with open(input_file) as fh:
        data = csv.DictReader(fh)
        for record in data:
            final_draft = {}

            final_draft["locus"] = record["gene symbol"]
            final_draft["panels"] = [panel_name]
            final_draft["public_comment"] = ""
            final_draft["private_comment"] = ""
            final_draft["cross_cutting_modifier"] = []
            final_draft["variant_types"] = []
            final_draft["variant_descriptions"] = []
            final_draft["variant_consequences"] = []
            final_draft["confidence"] = ""
            final_draft["phenotypes"] = []
            final_draft["allelic_requirement"] = record["allelic requirement"]

            # Publications
            final_draft["publications"] = []
            if record["g2p id"] in g2p_to_pmids:
                for pmid in g2p_to_pmids[record["g2p id"]]:
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

            # Mechanism
            final_draft["molecular_mechanism"] = {
                "name": "",
                "support": "",
            }
            final_draft["mechanism_synopsis"] = []
            # Mechanism evidence
            final_draft["mechanism_evidence"] = []

            # Disease
            final_draft["disease"] = {
                "disease_name": record["disease"],
                "cross_references": [],
            }

            # Build the session name using gene, allelic requirement and disease to check for duplicates before inserting the draft
            final_draft["session_name"] = build_session_name(
                record["gene symbol"],
                record["allelic requirement"],
                final_draft["disease"]["disease_name"],
            )

            # Check if the draft already exists based on the session name
            if final_draft["session_name"] in automatic_drafts:
                print(
                    f"Draft record for session '{final_draft['session_name']}' already exists. Skipping insertion."
                )
                continue

            final_draft["source_data"] = {"name": source}

            if not final_draft["publications"]:
                print(f"No publications found for {record['gene symbol']}.")
                continue

            # ### TO REMOVE ###
            # if record['gene symbol'] == "CRYAB":
            #     continue

            print(f"Creating draft for {record['gene symbol']}")

            # Call G2P API to insert the curation draft
            handle_draft_record(
                final_draft,
                api_url,
                cookies,
                output_records,
            )


def prepare_publications_pipeline(pipeline_file):
    """
    """
    g2p_to_pmids = {}

    with open(pipeline_file) as fh:
        data = csv.DictReader(fh)
        for record in data:
            pmid = record["PMID"]
            g2p_id = record["G2P_IDs"]

            if "gemini_relevance_label" in record:
                score = record["gemini_relevance_label"]

                if score != "high":
                    continue

            if g2p_id not in g2p_to_pmids:
                g2p_to_pmids[g2p_id] = []
            g2p_to_pmids[g2p_id].append(pmid)

    return g2p_to_pmids


def prepare_draft_records_immuno(
    input_file: str,
    immuno_file: str,
    panel_name: str,
    source: str,
    automatic_drafts: dict,
    api_url: str,
    cookies: requests.cookies.RequestsCookieJar,
    output_records: list | None = None,
) -> None:
    """
    Read the input file from the immuno panel to build the draft records.
    It also attaches ClinGen data to the draft record if the gene is present in
    the ClinGen input file.
    """
    import openpyxl

    ar_mapping = {
        "AR": "biallelic_autosomal",
        "AD": "monoallelic_autosomal",
        "XL": "monoallelic_X",
    }

    mechanism_mapping = {
        "GOF": "gain of function",
        "LOF": "loss of function",
        "DN": "dominant negative",
    }

    if not os.path.isfile(immuno_file):
        sys.exit(f"Invalid immuno file '{immuno_file}'")

    workbook = openpyxl.load_workbook(immuno_file, read_only=True, data_only=True)
    worksheet = workbook[workbook.sheetnames[0]]

    rows = worksheet.iter_rows(values_only=True)
    try:
        raw_headers = next(rows)
    except StopIteration:
        sys.exit(f"Immuno file '{immuno_file}' is empty")

    headers = normalize_excel_headers(raw_headers)
    immuno_records = []

    for row in rows:
        if row is None or not any(value is not None and str(value).strip() != "" for value in row):
            continue

        record = {}
        for index, header in enumerate(headers):
            value = row[index] if index < len(row) else None
            if isinstance(value, str):
                value = value.strip()
            record[header] = value

        # print("\n\n->", record)

        if "G2P record" in record and record["G2P record"] is not None:
            print(f"Gene is already in G2P: {record['Genetic defect']}. Skipping draft creation.")
            continue

        disease_name = record["Disease"].strip()
        gene_symbol = record["Genetic defect"].strip()

        if record["Inheritance"] is None:
            allelic_requirement = ""
        else:
            allelic_requirement = ar_mapping.get(record["Inheritance"].strip(), "")

        if record["GOF/DN"] is None:
            mechanism = ""
        else:
            mechanism = mechanism_mapping.get(record["GOF/DN"].strip(), "")

        # Check if mechanism is in disease name
        if mechanism == "":
            if "(LOF)" in disease_name or " LOF" in disease_name:
                mechanism = "loss of function"
            elif "(GOF)" in disease_name or " GOF" in disease_name:
                mechanism = "gain of function"
            elif "(DN)" in disease_name or " DN" in disease_name:
                mechanism = "dominant negative"

        # Attach data to record
        record["allelic_requirement"] = allelic_requirement
        record["mechanism"] = mechanism
        record["gene_symbol"] = gene_symbol
        record["disease_name"] = disease_name

        # print("---> Gene:", gene_symbol, "Disease:", disease_name, "AR:", allelic_requirement, "Mechanism:", mechanism)

        immuno_records.append(record)

    workbook.close()

    pre_draft_dict = {}
    with open(input_file) as fh:
        pre_draft_data = json.load(fh)
        for record in pre_draft_data:
            if record["gene_symbol"] not in pre_draft_dict:
                pre_draft_dict[record["gene_symbol"]] = [record]
            else:
                pre_draft_dict[record["gene_symbol"]].append(record)

    count_gene_symbols = {}

    for immuno_record in immuno_records:
        if immuno_record["gene_symbol"] not in pre_draft_dict:
            continue

        if immuno_record["gene_symbol"] in count_gene_symbols:
            count_gene_symbols[immuno_record["gene_symbol"]] += 1
        else:
            count_gene_symbols[immuno_record["gene_symbol"]] = 1

        disease_name = immuno_record["disease_name"]
        mechanism = immuno_record["mechanism"]
        pre_drafts_clingen = pre_draft_dict[immuno_record["gene_symbol"]]
        clingen_record = None

        if len(pre_drafts_clingen) > 1:
            for pre_draft in pre_drafts_clingen:
                ar = allelic_requirement_mapping.get(pre_draft["allelic_requirement"], "")
                if ar == immuno_record["allelic_requirement"] and (
                    pre_draft["mechanism"] == mechanism
                    or pre_draft["mechanism"] == ""
                    or mechanism == ""
                ):
                    clingen_record = pre_draft
                    break
        else:
            if count_gene_symbols[immuno_record["gene_symbol"]] > 1:
                print(f"WARNING: Multiple immuno rows for gene: {immuno_record['gene_symbol']}")

            ar = allelic_requirement_mapping.get(pre_drafts_clingen[0]["allelic_requirement"], "")
            if ar == immuno_record["allelic_requirement"] and (
                pre_drafts_clingen[0]["mechanism"] == mechanism
                or pre_drafts_clingen[0]["mechanism"] == ""
                or mechanism == ""
            ):
                clingen_record = pre_drafts_clingen[0]
            else:
                print(f"WARNING: No matching allelic requirement/mechanism for gene: {immuno_record['gene_symbol']}. Skipping.")
                continue

        final_draft = {}
        final_draft["locus"] = immuno_record["gene_symbol"]
        final_draft["panels"] = [panel_name]
        final_draft["public_comment"] = ""
        final_draft["private_comment"] = ""
        final_draft["cross_cutting_modifier"] = []
        final_draft["variant_types"] = []
        final_draft["variant_descriptions"] = []
        final_draft["variant_consequences"] = []
        final_draft["confidence"] = ""
        final_draft["phenotypes"] = []
        final_draft["allelic_requirement"] = immuno_record["allelic_requirement"]

        if clingen_record is None:
            print(f"WARNING: No matching ClinGen record found for gene: {immuno_record['gene_symbol']}, allelic requirement: {immuno_record['allelic_requirement']}, mechanism: {mechanism}. Skipping.")
            continue

        final_draft["publications"] = []
        for pmid in clingen_record["pmids"]:
            publication_data = fetch_pmid_info(api_url, pmid)
            if publication_data["authors"] is not None:
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

        final_draft["molecular_mechanism"] = {
            "name": mechanism,
            "support": "",
        }
        final_draft["mechanism_synopsis"] = []
        final_draft["mechanism_evidence"] = []

        final_draft["disease"] = {
            "disease_name": immuno_record["disease_name"],
            "cross_references": [],
        }

        final_draft["source_data"] = {"name": source}
        if "url" in clingen_record:
            final_draft["source_data"]["url"] = clingen_record["url"]

        final_draft["source_data"]["mechanism"] = clingen_record["mechanism"]
        mechanism_evidence = build_clingen_mechanism_evidence(clingen_record)
        if mechanism_evidence:
            final_draft["source_data"]["mechanism_evidence"] = mechanism_evidence

        final_draft["source_data"]["phenotypes"] = clingen_record["phenotypes"]

        disease_cross_references = []
        if "mondo_id" in clingen_record and clingen_record["mondo_id"] != "":
            disease_data = fetch_disease_info(api_url, clingen_record["mondo_id"])
            disease_cross_references.append(
                {
                    "source": "Mondo",
                    "identifier": clingen_record["mondo_id"],
                    "disease_name": disease_data["disease"].lower(),
                    "original_disease_name": disease_data["disease"],
                }
            )
        if "disease_id" in clingen_record and clingen_record["disease_id"] != "" and clingen_record["disease_id"] != []:
            omim_list = clingen_record["disease_id"].split(",")
            for omim_id in omim_list:
                if omim_id.startswith("OMIM:") or omim_id.startswith("MIM:"):
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
                                "original_disease_name": disease_data["disease"],
                            }
                        )
                elif omim_id.startswith("MONDO:"):
                    mondo_id = omim_id.strip()
                    disease_data = fetch_disease_info(api_url, mondo_id)
                    disease_cross_references.append(
                        {
                            "source": "Mondo",
                            "identifier": mondo_id,
                            "disease_name": disease_data["disease"].lower(),
                            "original_disease_name": disease_data["disease"],
                        }
                    )

        # Save the disease name and cross references in the source field
        if "disease" in clingen_record:
            final_draft["source_data"]["disease"] = clingen_record["disease"]
        final_draft["source_data"]["disease_cross_references"] = (
            disease_cross_references
        )

        final_draft["session_name"] = build_session_name(
            immuno_record["gene_symbol"],
            immuno_record["allelic_requirement"],
            immuno_record["disease_name"],
        )

        # Check if the draft already exists based on the session name
        if final_draft["session_name"] in automatic_drafts:
            print(
                f"Draft record for session '{final_draft['session_name']}' already exists. Skipping insertion."
            )
            continue

        # print(f"Creating draft for {immuno_record['gene_symbol']} with disease {immuno_record['disease_name']}")

        # print(json.dumps(final_draft, indent=4))

        # Call G2P API to insert the curation draft
        handle_draft_record(
            final_draft,
            api_url,
            cookies,
            output_records,
        )


def main():
    parser = argparse.ArgumentParser(description="Create draft records")
    parser.add_argument(
        "--source",
        required=True,
        choices=["clingen", "immuno", "inhouse"],
        help="Draft source: clingen, immuno, or inhouse",
    )
    parser.add_argument(
        "--input_file",
        required=True,
        help="Source-specific input file: ClinGen JSON, or inhouse CSV. For immuno this is the ClinGen JSON.",
    )
    parser.add_argument(
        "--immuno_file",
        required=False,
        help="Required when --source is immuno. Immunology workbook (xlsx).",
    )
    parser.add_argument(
        "--clingen_panel",
        required=False,
        help="ClinGen panel to process when --source is clingen",
    )
    parser.add_argument(
        "--genes_to_include",
        required=False,
        help="File containing gene symbols to include, one per line",
    )
    parser.add_argument(
        "--config", required=True, help="Config file with details to G2P API"
    )
    parser.add_argument(
        "--panel", required=True, help="G2P panel name e.g. 'Ear disorders'"
    )
    parser.add_argument(
        "--output_file", required=False, help="Output file to save the final draft records"
    )
    args = parser.parse_args()

    input_file = args.input_file
    source = args.source
    panel_name = args.panel
    output_file = args.output_file
    immuno_file = args.immuno_file
    clingen_panel = args.clingen_panel
    genes_to_include = args.genes_to_include
    normalized_source = source.lower()

    if normalized_source == "immuno" and not immuno_file:
        sys.exit("ERROR: --immuno_file is required when --source is immuno")
    if normalized_source != "immuno" and immuno_file:
        sys.exit("ERROR: --immuno_file can only be used when --source is immuno")

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
    output_records = [] if output_file else None

    if normalized_source == "inhouse":
        # Prepare the publications
        # Read the output file generated by the pipeline
        print(f"Preparing drafts from inhouse file: {input_file}...")
        g2p_to_pmids = prepare_publications_pipeline(input_file)

        prepare_draft_records_pipeline(
            input_file, g2p_to_pmids, panel_name, normalized_source, automatic_drafts, api_url, cookies, output_records
        )
        print(f"Preparing drafts from inhouse file: {input_file}... done")
    elif normalized_source == "immuno":
        # Prepare the immuno panel genes
        print(f"Preparing drafts from immuno file: {immuno_file}...")
        prepare_draft_records_immuno(
            input_file, immuno_file, panel_name, normalized_source, automatic_drafts, api_url, cookies, output_records
        )
        print(f"Preparing drafts from immuno file: {immuno_file}... done")
    else:
        clingen_panel_message = (
            f"ClinGen panel: {clingen_panel}" if clingen_panel else "all ClinGen panels"
        )
        print(
            f"Preparing ClinGen drafts from Gemini input file: {input_file} "
            f"and {clingen_panel_message}..."
        )
        prepare_clingen_draft_records(
            input_file,
            panel_name,
            normalized_source,
            automatic_drafts,
            api_url,
            cookies,
            clingen_panel,
            genes_to_include,
            output_records,
        )
        print(
            f"Preparing ClinGen drafts from Gemini input file: {input_file} "
            f"and {clingen_panel_message}... done"
        )

    if output_file:
        write_draft_output(output_file, output_records)
        print(f"Draft records written to {output_file}")

    logout(api_url, cookies)


if __name__ == "__main__":
    main()
