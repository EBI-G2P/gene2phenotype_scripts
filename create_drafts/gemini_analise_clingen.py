#!/usr/bin/env python3

import argparse
import configparser
import json
import os
import sys
from pathlib import Path
from google import genai
from google.genai.types import HttpOptions
from google.oauth2 import service_account
from pydantic import BaseModel


"""
Script to extract data from the ClinGen evidence summaries using the Google Gemini model.

Options:
        --config:        Config file containing the Gemini model details (mandatory)
                            File format is the following: 
                            [project_config]
                            key_file = <> # Google Vertex AI key (format: json)
                            model = <>    # Gemini model (default: gemini-2.5-flash)
                            location = <> # Location for the Gemini model (default: europe-west2)
                            project = <>
        --input_file:    Input JSON file with ClinGen evidence summaries (mandatory)
        --output_file:   Output JSON file with Gemini analysis (default: <input_stem>_gemini.json)
        -l/--limit:      Process N entries and exit (default: all)
        --clingen_panel: Run the analysis only for records in the specific ClinGen panel (default: all)

The input file is read-only. The script writes analysed records to the output file.
If the output file already exists, records already present in that file are skipped
so interrupted runs can be resumed.

Example usage:
    python gemini_analise_clingen.py \
        --input_file clingen_extracted_data_2026-09-25.json \
        --output_file clingen_extracted_data_2026-09-25_gemini.json \
        --config config.ini \
        --clingen_panel "Hearing Loss Gene Curation Expert Panel"
"""


def run_process(args):
    with args.input_file.open("rt") as fh:
        clingen_data = json.load(fh)

    output_file = args.output_file or default_output_file(args.input_file)
    analysed_data = load_existing_output(output_file)
    analysed_record_keys = {record_key(record) for record in analysed_data}
    records_to_process = [
        record
        for record in clingen_data
        if record_key(record) not in analysed_record_keys
        and (not args.clingen_panel or args.clingen_panel == record["clingen_panel"])
    ]

    if not records_to_process:
        write_output(output_file, analysed_data)
        return

    # Get the Gemini model details from the config file
    config = configparser.ConfigParser()
    config.read(args.config)
    try:
        gemini_config = config["project_config"]
    except KeyError:
        sys.exit("ERROR: 'project_config' missing from config file")
    args.key_file = gemini_config.get("key_file")
    args.model = gemini_config.get("model", "gemini-2.5-flash")
    args.location = gemini_config.get("location", "europe-west2")
    args.project_name = gemini_config.get("project")

    credentials = load_json_key(args.key_file)

    client = genai.Client(
        vertexai=True,
        project=args.project_name,
        location=args.location,  # gemini pro is available at us-central1; flash is in europe-west2
        credentials=credentials,
        http_options=HttpOptions(api_version="v1"),
    )

    try:
        done = 0

        for record in records_to_process:
            output = process_publication(client, record, args.model)
            record["pmids"] = output.pmids
            record["disease_id"] = output.disease_id
            record["mechanism"] = output.mechanism
            record["allelic_requirement"] = output.allelic_requirement
            record["phenotypes"] = output.phenotypes
            record["evidence"] = output.experimental_evidence
            record["comment"] = output.comment
            analysed_data.append(record)
            analysed_record_keys.add(record_key(record))

            print(
                f"\nClinGen record gene: {record['gene_symbol']}, disease: {record['disease']}",
                file=sys.stderr,
            )
            print(f"Publications        : {output.pmids}", file=sys.stderr)
            print(f"Gene                : {output.gene}", file=sys.stderr)
            print(f"Disease             : {output.disease}", file=sys.stderr)
            print(f"OMIM/Mondo ID       : {output.disease_id}", file=sys.stderr)
            print(f"Mechanism           : {output.mechanism}", file=sys.stderr)
            print(
                f"Allelic requirement : {output.allelic_requirement}", file=sys.stderr
            )
            print(f"Phenotypes          : {output.phenotypes}", file=sys.stderr)
            print(f"Comment             : {output.comment}", file=sys.stderr)

            done += 1
            write_output(output_file, analysed_data)
            if done == args.limit:
                return
    finally:
        write_output(output_file, analysed_data)


def default_output_file(input_file: Path) -> Path:
    return input_file.with_name(f"{input_file.stem}_gemini{input_file.suffix}")


def load_existing_output(output_file: Path) -> list:
    if not output_file.is_file():
        return []

    with output_file.open("rt") as fh:
        return json.load(fh)


def write_output(output_file: Path, data: list) -> None:
    temp_output_file = output_file.with_name(f"{output_file.name}.tmp")
    with temp_output_file.open("wt") as fh:
        json.dump(data, fh, indent=2)
    os.replace(temp_output_file, output_file)


def record_key(record: dict) -> tuple:
    return (
        record.get("gene_symbol", ""),
        record.get("disease", ""),
        record.get("mondo_id", ""),
        record.get("clingen_panel", ""),
        record.get("url", ""),
    )


def load_json_key(key_file):
    credentials = service_account.Credentials.from_service_account_file(
        key_file
    ).with_scopes(["https://www.googleapis.com/auth/cloud-platform"])
    return credentials


class Relevance(BaseModel):
    pmids: list
    disease: str
    disease_id: str
    mechanism: str
    allelic_requirement: str
    gene: str
    phenotypes: list
    experimental_evidence: list
    comment: str


def process_publication(client: genai.Client, record: dict, model: str) -> Relevance:
    prompt = f"""\
You are a biomedical information extraction assistant. \
You will be provided with a specific gene and disease, along with an Evidence summary \
that contains scientific evidence for this gene-disease association.

Your task is to extract structured information only for this gene-disease pair.\
For this association extract: \
- pmids: A list of all PubMed IDs associated with the specific gene-disease. \
- disease: The disease or diseases mentioned. \
- disease_id: The MIM/OMIM or Mondo IDs associated with the \
specific disease. \
- mechanism: The mechanism of the specific disease if mentioned \
(e.g. "gain of function", "loss of function") \
- allelic_requirement: The allelic requirement or mode \
of inheritance if mentioned (e.g., "autosomal dominant", "autosomal recessive"). \
- gene: The gene symbol mentioned. \
- phenotypes: Any specific phenotypes described in the text
associated with the specific disease. \
(e.g. "reduced tendon reflexes", "distal motor weakness", "sensory disturbances"). \
- experimental_evidence: Any experimental evidence supporting \
the rule of the specific gene in the specific disease. The type of evidence \
can be: function (evidence related to gene expression or computer simulations), \
rescue (evidence showing that the phenotype can be rescued), models (a model
with a disrupted copy of the gene shows a phenotype consistent with the human disease) \
or functional alteration (evidence showing that cultured cells, in which the \
function of the gene has been disrupted, have a phenotype that is consistent \
with the human disease process). 

Only extract information that is explicitly mentioned in the text. \
If a field is not mentioned, return an empty list for pmids or phenotypes, \
and an empty string for gene, disease, disease_id, mechanism or allelic_requirement.
The specific gene-disease are provided in the input as gene and disease.
If there are PMIDs associated with other diseases do not include them \
in the output.

After producing the structured output, provide one short comment stating \
whether additional diseases (besides the input disease) appear in the evidence summary.

Input:
gene: {record["gene_symbol"]}\
disease: {record["disease"]}
Here it is the text to analise:
Evidence summary: {record["evidence_summary"]}\
"""
    response = client.models.generate_content(
        model=model,
        contents=prompt,
        config={
            "response_mime_type": "application/json",
            "response_schema": Relevance,
            "temperature": 0,
        },
    )

    return response.parsed


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--input_file",
        required=True,
        type=Path,
        help="Input JSON file with ClinGen evidence summaries"
    )
    parser.add_argument(
        "--output_file",
        type=Path,
        default=None,
        help="Output JSON file with Gemini analysis (default: <input_stem>_gemini.json)",
    )
    parser.add_argument(
        "--config",
        type=Path,
        required=True,
        help="Config file containing the Gemini model details (key_file, model, location, project_name)",
    )
    parser.add_argument(
        "-l",
        "--limit",
        type=int,
        default=0,
        metavar="N",
        help="Process N entries and exit (default: all)",
    )
    parser.add_argument(
        "--clingen_panel",
        type=str,
        default=None,
        help="Run the analysis only for records in the specific ClinGen panel",
    )

    args = parser.parse_args()

    run_process(args)


if __name__ == "__main__":
    main()
