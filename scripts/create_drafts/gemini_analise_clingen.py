#!/usr/bin/env python3

import argparse
import json
import sys
import time
from pathlib import Path
from google import genai
from google.genai.types import HttpOptions
from google.oauth2 import service_account
from pydantic import BaseModel


"""

python scripts/create_drafts/gemini_analise_clingen.py \
    --input_file clingen_extracted_data.json \
    --key_file project_key.json \
    -m gemini-2.5-pro --location us-central1 \
    --clingen_panel "Hearing Loss Gene Curation Expert Panel"
"""


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--input_file", type=Path, help="Input JSON file")
    parser.add_argument(
        "--key_file",
        type=Path,
        required=True,
        help="Google Vertex AI key (json file)"
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
        "-m",
        "--model",
        default="gemini-2.5-flash",
        help="Gemini model (default: gemini-2.5-flash)",
    )
    parser.add_argument(
        "--location",
        default="europe-west2",
        help="Location for the Gemini model (default: europe-west2)",
    )
    parser.add_argument(
        "--rpm",
        type=int,
        default=10,
        help="Max requests per minute (default: 10)"
    )
    parser.add_argument(
        "--clingen_panel",
        type=str,
        default=None,
        help="Run the analysis only for records in the specific ClinGen panel"
    )

    args = parser.parse_args()

    run_process(args)

def run_process(args):
    with args.input_file.open("rt") as fh:
        clingen_data = json.load(fh)

    credentials = load_json_key(args.key_file)

    client = genai.Client(
        vertexai=True,
        project="prj-int-dev-paradigm-g2p",
        location=args.location, # gemini pro is available at us-central1; flash is in europe-west2
        credentials=credentials,
        http_options=HttpOptions(api_version="v1")
    )

    try:
        done = 0
        
        if args.rpm > 0:
            secs = 60.0 / args.rpm
        else:
            secs = 0

        for record in clingen_data:
            if "pmids" in record:
                continue

            if args.clingen_panel:
                if args.clingen_panel != record["clingen_panel"]:
                    continue

            output = process_publication(client, record, args.model)
            record["pmids"] = output.pmids
            record["disease_id"] = output.disease_id
            record["mechanism"] = output.mechanism
            record["allelic_requirement"] = output.allelic_requirement
            record["phenotypes"] = output.phenotypes
            record["evidence"] = output.experimental_evidence
            record["comment"] = output.comment

            print(
                f"\nClinGen record gene: {record['gene_symbol']}, disease: {record['disease']}",
                file=sys.stderr,
            )
            print(f"Publications        : {output.pmids}", file=sys.stderr)
            print(f"Gene                : {output.gene}", file=sys.stderr)
            print(f"Disease             : {output.disease}", file=sys.stderr)
            print(f"OMIM/Mondo ID       : {output.disease_id}", file=sys.stderr)
            print(f"Mechanism           : {output.mechanism}", file=sys.stderr)
            print(f"Allelic requirement : {output.allelic_requirement}", file=sys.stderr)
            print(f"Phenotypes          : {output.phenotypes}", file=sys.stderr)
            print(f"Comment             : {output.comment}", file=sys.stderr)

            done += 1
            if done == args.limit:
                return
            
            time.sleep(secs)
    finally:
        with args.input_file.open("wt") as fh:
            json.dump(clingen_data, fh, indent=2)


def load_json_key(key_file):
    credentials = service_account.Credentials.from_service_account_file(key_file).with_scopes(["https://www.googleapis.com/auth/cloud-platform"])
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


def process_publication(
    client: genai.Client, record: dict, model: str
) -> Relevance:
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
gene: {record['gene_symbol']}\
disease: {record['disease']}
Here it is the text to analise:
Evidence summary: {record['evidence_summary']}\
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


if __name__ == "__main__":
    main()
