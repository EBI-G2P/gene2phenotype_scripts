#!/usr/bin/env python3

import sys
import argparse
import configparser
import csv
import os.path
from html.parser import HTMLParser
import requests
import json
from typing import List, Dict


"""
Script to extract the ClinGen evidence summary.

Options:
        --input_file:    ClinGen gene-disease validity (supported format: csv)

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
        if text == "Evidence Summary:":
            self.just_saw_label = True
            return

        # Collect text inside summary cell
        if self.in_summary_cell:
            self.current_summary.append(text)


def load_data(file: str) -> List[Dict[str, str]]:
    """
    Read the ClinGen file
    """
    with open(file, "r", encoding='utf-8') as fh:
        header = []

        for line in fh:
            if line.startswith('"GENE SYMBOL"') or line.startswith("GENE SYMBOL"):
                header = line
                break      

        # Create DictReader using the header
        reader = csv.DictReader([header] + list(fh), skipinitialspace=True)

        # Read rows skipping separator rows (++++++++)
        rows = []
        for row in reader:
            if row["GENE SYMBOL"].startswith("+"):
                continue
            rows.append(row)

        return rows


def process_clingen_data(clingen_data: List[Dict[str, str]]) -> None:
    output_file = "clingen_extracted_data.json"
    final_data = []

    for row in clingen_data:
        url = row["ONLINE REPORT"]

        response = requests.get(url)
        parser = EvidenceSummaryParser()
        parser.feed(response.text)
        evidence_summary = parser.all_summaries

        final_data.append({
            "gene_symbol": row["GENE SYMBOL"],
            "hgnc_id": row["GENE ID (HGNC)"],
            "disease": row["DISEASE LABEL"],
            "mondo_id": row["DISEASE ID (MONDO)"],
            "confidence": row["CLASSIFICATION"],
            "clingen_panel": row["GCEP"],
            "evidence_summary": evidence_summary,
            "url": url
        })

    with open(output_file, "w", encoding="utf-8") as fw:
        json.dump(final_data, fw, ensure_ascii=False, indent=2)


def main():
    parser = argparse.ArgumentParser(
        description="Fetch gene-disease data from ClinGen"
    )
    parser.add_argument(
        "--input_file",
        required=True,
        help="ClinGen gene-disease validity input file (csv)"
    )
    args = parser.parse_args()

    input_file = args.input_file

    if not os.path.isfile(input_file):
        sys.exit(f"Invalid input file '{input_file}'")

    clingen_data = load_data(input_file)

    # Extract the evidence summary from the ClinGen URL
    # Write output data to a json file
    process_clingen_data(clingen_data)

if __name__ == "__main__":
    main()
