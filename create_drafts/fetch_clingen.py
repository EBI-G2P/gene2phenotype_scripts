#!/usr/bin/env python3

import argparse
import csv
import datetime
import json
import os.path
import sys
from typing import List, Dict, Optional, Tuple
from html.parser import HTMLParser
import requests

REQUEST_TIMEOUT_SECONDS = 30
DEFAULT_OUTPUT_FILE_PREFIX = "clingen_extracted_data"
REQUIRED_INPUT_COLUMNS = {
    "GENE SYMBOL",
    "GENE ID (HGNC)",
    "DISEASE LABEL",
    "DISEASE ID (MONDO)",
    "CLASSIFICATION",
    "GCEP",
    "ONLINE REPORT",
}
EXPERIMENTAL_EVIDENCE_FIELD_ALIASES = {
    "label": "label",
    "experimental category": "experimental_category",
    "reference": "reference",
    "explanation": "explanation",
    "reason for changed score": "reason_for_changed_score",
    "evidence type": "evidence_type",
}


"""
Script to extract the ClinGen evidence summary and experimental evidence.

Options:
        --input_file (mandatory):   ClinGen gene-disease validity (supported format: csv)
                                    Download from: https://search.clinicalgenome.org/kb/downloads#section_gene-disease-validity
                                    File section: "Gene-Disease Validity Curations"

        --output_file (optional):   Output JSON file (default: clingen_extracted_data_YYYY-MM-DD.json)

Example how to run:
    python fetch_clingen.py --input_file Clingen-Gene-Disease-Summary-2025-11-11.csv
"""


class EvidenceSummaryParser(HTMLParser):
    """
    Class that extends HTMLParser to fetch evidence summary
    from the ClinGen HTML page.
    """
    def __init__(self):
        super().__init__()
        self.just_saw_label = False
        self.in_summary_cell = False
        self.current_summary = []
        self.all_summaries = []

    def handle_starttag(self, tag, attrs):
        if self.just_saw_label and tag == "td":
            # Start of the summary cell
            self.in_summary_cell = True
            # Reset the saw label
            self.just_saw_label = False

    def handle_endtag(self, tag):
        if tag == "td" and self.in_summary_cell:
            # Finished one summary cell
            self.in_summary_cell = False
            summary_text = " ".join(self.current_summary).strip()
            if summary_text:
                self.all_summaries.append(summary_text)
            self.current_summary = []

    def handle_data(self, data):
        text = data.strip()
        if not text:
            return

        # Detect the label containing the evidence we want to extract
        if _normalize_cell_value(text).casefold().rstrip(":").strip() == "evidence summary":
            self.just_saw_label = True
            return

        # Collect text inside summary cell
        if self.in_summary_cell:
            self.current_summary.append(text)


class ExperimentalEvidenceTableParser(HTMLParser):
    """
    Class that extends HTMLParser to extract the Experimental Evidence
    table from the ClinGen HTML page.
    """
    def __init__(self):
        super().__init__()
        self.in_heading = False
        self.heading_parts = []
        self.in_section = False
        self.in_table = False
        self.tables = []
        self.current_table = []
        self.in_row = False
        self.current_row = []
        self.in_cell = False
        self.current_cell = []

    def handle_starttag(self, tag, attrs):
        if tag in {"h1", "h2", "h3", "h4"}:
            self.in_heading = True
            self.heading_parts = []

        if self.in_section and tag == "table":
            self.in_table = True
            self.current_table = []

        if self.in_table and tag == "tr":
            self.in_row = True
            self.current_row = []

        if self.in_row and tag in {"td", "th"}:
            self.in_cell = True
            self.current_cell = []

        if self.in_cell and tag == "br":
            self.current_cell.append("\n")

    def handle_endtag(self, tag):
        if tag in {"h1", "h2", "h3", "h4"} and self.in_heading:
            heading_text = " ".join(self.heading_parts).strip()
            normalized_heading = _normalize_cell_value(heading_text).casefold()

            if normalized_heading.startswith("experimental evidence"):
                self.in_section = True
            elif self.in_section and normalized_heading.startswith("total points"):
                self.in_section = False

            self.in_heading = False
            self.heading_parts = []

        if tag in {"td", "th"} and self.in_cell:
            cell_text = " ".join(" ".join(self.current_cell).split())
            self.current_row.append(cell_text)
            self.current_cell = []
            self.in_cell = False

        if tag == "tr" and self.in_row:
            if self.current_row:
                self.current_table.append(self.current_row)
            self.current_row = []
            self.in_row = False

        if tag == "table" and self.in_table:
            if self.current_table:
                self.tables.append(self.current_table)
            self.current_table = []
            self.in_table = False

    def handle_data(self, data):
        text = data.strip()
        if not text:
            return

        if self.in_heading:
            self.heading_parts.append(text)
            return

        if self.in_cell:
            self.current_cell.append(text)


def load_data(file: str) -> List[Dict[str, str]]:
    """
    Read the ClinGen file
    """
    with open(file, "r", encoding='utf-8') as fh:
        header = ""

        for line in fh:
            if line.startswith('"GENE SYMBOL"') or line.startswith("GENE SYMBOL"):
                header = line
                break      

        if not header:
            sys.exit("ERROR: Could not find ClinGen CSV header row starting with 'GENE SYMBOL'")

        # Create DictReader using the header
        reader = csv.DictReader([header] + list(fh), skipinitialspace=True)
        missing_columns = REQUIRED_INPUT_COLUMNS - set(reader.fieldnames or [])
        if missing_columns:
            missing_columns_list = ", ".join(sorted(missing_columns))
            sys.exit(f"ERROR: ClinGen CSV is missing required columns: {missing_columns_list}")

        # Read rows skipping separator rows (++++++++)
        rows = []
        for row in reader:
            if row["GENE SYMBOL"].startswith("+"):
                continue
            rows.append(row)

        return rows


def _normalize_cell_value(text: str) -> str:
    return " ".join(text.split()).strip()


def _normalize_header(text: str) -> str:
    return _normalize_cell_value(text).lower()


def _find_experimental_evidence_table(
    tables: List[List[List[str]]]
) -> Tuple[Optional[List[List[str]]], Optional[int]]:
    for table in tables:
        for row_index, row in enumerate(table):
            normalized_row = [_normalize_header(cell) for cell in row]
            if (
                "label" in normalized_row
                and "experimental category" in normalized_row
                and "reference" in normalized_row
            ):
                return table, row_index

    return None, None


def extract_experimental_evidence(html_text: str) -> List[Dict[str, str]]:
    parser = ExperimentalEvidenceTableParser()
    parser.feed(html_text)

    table, header_row_index = _find_experimental_evidence_table(parser.tables)
    if table is None or header_row_index is None:
        return []

    header_row = table[header_row_index]
    column_indexes = {}
    for idx, cell in enumerate(header_row):
        normalized_header = _normalize_header(cell)
        field_name = EXPERIMENTAL_EVIDENCE_FIELD_ALIASES.get(normalized_header)
        if field_name:
            column_indexes[field_name] = idx

    required_fields = {
        "label",
        "experimental_category",
        "reference",
        "explanation",
        "reason_for_changed_score",
    }
    if not required_fields.issubset(column_indexes):
        return []

    records = []
    for row in table[header_row_index + 1:]:
        if not any(_normalize_cell_value(cell) for cell in row):
            continue

        first_cell = _normalize_header(row[0])
        if first_cell.startswith("total points"):
            break

        record = {}
        for field_name, idx in column_indexes.items():
            value = row[idx] if idx < len(row) else ""
            record[field_name] = _normalize_cell_value(value)

        if not record.get("label"):
            continue

        records.append(record)

    return records


def default_output_file() -> str:
    today = datetime.date.today().isoformat()
    return f"{DEFAULT_OUTPUT_FILE_PREFIX}_{today}.json"


def write_output(output_file: str, final_data: List[Dict[str, str]]) -> None:
    temp_output_file = f"{output_file}.tmp"
    with open(temp_output_file, "w", encoding="utf-8") as fw:
        json.dump(final_data, fw, ensure_ascii=False, indent=2)
    os.replace(temp_output_file, output_file)


def warn_empty_data(
    gene_symbol: str,
    disease: str,
    url: str,
    evidence_summary: List[str],
    experimental_evidence: List[Dict[str, str]],
) -> None:
    context = f"{gene_symbol} / {disease} ({url})"
    if not evidence_summary:
        print(f"WARNING: No evidence summary extracted for {context}", file=sys.stderr)
    if not experimental_evidence:
        print(f"WARNING: No experimental evidence extracted for {context}", file=sys.stderr)


def process_clingen_data(clingen_data: List[Dict[str, str]], output_file: str) -> None:
    final_data = []

    for row in clingen_data:
        url = row["ONLINE REPORT"]
        gene_symbol = row["GENE SYMBOL"]
        disease = row["DISEASE LABEL"]

        try:
            response = requests.get(url, timeout=REQUEST_TIMEOUT_SECONDS)
            response.raise_for_status()
        except requests.exceptions.RequestException as error:
            print(
                f"WARNING: Skipping {gene_symbol} / {disease} ({url}): {error}",
                file=sys.stderr,
            )
            write_output(output_file, final_data)
            continue

        parser = EvidenceSummaryParser()
        parser.feed(response.text)
        evidence_summary = parser.all_summaries
        experimental_evidence = extract_experimental_evidence(response.text)
        warn_empty_data(gene_symbol, disease, url, evidence_summary, experimental_evidence)

        final_data.append({
            "gene_symbol": gene_symbol,
            "hgnc_id": row["GENE ID (HGNC)"],
            "disease": disease,
            "mondo_id": row["DISEASE ID (MONDO)"],
            "confidence": row["CLASSIFICATION"],
            "clingen_panel": row["GCEP"],
            "evidence_summary": evidence_summary,
            "experimental_evidence": experimental_evidence,
            "url": url
        })
        write_output(output_file, final_data)

    write_output(output_file, final_data)


def main():
    parser = argparse.ArgumentParser(
        description="Fetch gene-disease data from ClinGen"
    )
    parser.add_argument(
        "--input_file",
        required=True,
        help="ClinGen gene-disease validity input file (csv)"
    )
    parser.add_argument(
        "--output_file",
        default=default_output_file(),
        help="Output JSON file (default: clingen_extracted_data_YYYY-MM-DD.json)"
    )
    args = parser.parse_args()

    input_file = args.input_file

    if not os.path.isfile(input_file):
        sys.exit(f"Invalid input file '{input_file}'")

    clingen_data = load_data(input_file)

    # Extract the evidence summary from the ClinGen URL
    # Write output data to a json file
    process_clingen_data(clingen_data, args.output_file)

if __name__ == "__main__":
    main()
