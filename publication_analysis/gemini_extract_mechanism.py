#!/usr/bin/env python3

import argparse
import csv
import configparser
import json
import sys
import time
from html.parser import HTMLParser
from io import StringIO
from pathlib import Path
from urllib.request import urlopen
from urllib.error import HTTPError
import xml.etree.ElementTree as ET

from google import genai
from google.genai.types import HttpOptions
from google.oauth2 import service_account
from pydantic import BaseModel


"""
Script to download G2P records and associated publications (mined and curated) and to
extract mechanism information for each publication using Google Vertex AI Gemini models.

Usage:
    1) Download G2P records and associated publications
    $ python gemini_extract_mechanism.py init g2p_records.json

    2) Extract mechanism for 300 publications
    $ python gemini_extract_mechanism.py process --config config.ini --limit 300 g2p_records.json

The config.ini file should contain the following entries:
    [project_config]
    key_file = path/to/key_file.json
    project = your-project-id
    location = europe-west2
    model = gemini-2.5-flash

Alternatively, you can provide a text file with G2P record IDs to analyse only those records:
    --file_records g2p_record_ids.txt
"""


def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(required=True)

    p1 = subparsers.add_parser("init")
    p1.add_argument("output", type=Path, help="Output JSON file")
    p1.set_defaults(func=run_download)
    # TODO: add another mode to append new records with mined publications
    p2 = subparsers.add_parser("process")
    p2.add_argument("infile", type=Path, help="Input JSON file")
    p2.add_argument(
        "--config",
        type=Path,
        required=True,
        help="Config file with Google Vertex AI settings(key file, project name, model, etc.)",
    )
    p2.add_argument(
        "--file_records",
        type=Path,
        required=False,
        help="G2P records to analyse (txt file)"
    )
    p2.add_argument(
        "-l",
        "--limit",
        type=int,
        default=0,
        metavar="N",
        help="Process N publications and exit (default: all)",
    )
    p2.add_argument(
        "--rpm",
        type=int,
        default=10,
        help="Max requests per minute (default: 10)"
    )
    p2.set_defaults(func=run_process)

    args = parser.parse_args()
    args.func(args)


def run_download(args):
    """
    Calls method to download the G2P records.
    Writes the records to a json file.
    """
    if args.output.is_file():
       sys.exit(f"File '{args.output}' already exists")

    with args.output.open("wt") as fh:
        records = download_g2p()
        json.dump(records, fh, indent=2)


def run_process(args):
    with args.infile.open("rt") as fh:
        records = json.load(fh)

    # Read config file
    config_parser = configparser.ConfigParser()
    config_parser.read(args.config)
    config = config_parser["project_config"]

    if "key_file" not in config or "project" not in config:
        sys.exit("Error: 'key_file' or 'project' not found in config file")

    credentials = load_json_key(config["key_file"])

    if "location" not in config:
        config["location"] = "europe-west2" # default location

    if "model" not in config:
        config["model"] = "gemini-2.5-flash" # default model

    client = genai.Client(
        vertexai=True,
        project=config["project"],
        location=config["location"], # gemini pro is available at us-central1; flash is in europe-west2
        credentials=credentials,
        http_options=HttpOptions(api_version="v1")
    )

    if args.file_records:
        with args.file_records.open("rt") as fh:
            records_to_analyse = set(line.strip() for line in fh if line.strip())

    try:
        done = 0
        
        if args.rpm > 0:
            secs = 60.0 / args.rpm
        else:
            secs = 0

        for record in records:
            for pub in record["publications"]:
                if records_to_analyse and record["id"] not in records_to_analyse:
                    continue

                if pub["mechanism"] is not None:
                    continue

                article = get_article(pub["id"])
                if article is None:
                    pub["mechanism"] = "incomplete"
                    continue

                pub.update(**article)

                output = process_publication(client, record, pub, config["model"])
                pub["mechanism"] = output.mechanism
                pub["mechanism_evidence"] = output.mechanism_evidence
                pub["comment"] = output.comment

                if done:
                    print("", file=sys.stderr)

                print(
                    f"G2P        : {record['id']}, {record['gene']}, {record['disease']}, {record['confidence']}",
                    file=sys.stderr,
                )
                print(f"Article  : {pub['title']}", file=sys.stderr)
                print(f"Mechanism: {pub['mechanism']}", file=sys.stderr)
                print(f"Evidence : {pub['mechanism_evidence']}", file=sys.stderr)
                print(f"Comment  : {pub['comment']}", file=sys.stderr)

                done += 1
                if done == args.limit:
                    return

                time.sleep(secs)
    finally:
        with args.infile.open("wt") as fh:
            json.dump(records, fh, indent=2)


def load_json_key(key_file):
    credentials = service_account.Credentials.from_service_account_file(key_file).with_scopes(["https://www.googleapis.com/auth/cloud-platform"])
    return credentials


class Relevance(BaseModel):
    mechanism: str
    mechanism_evidence: list
    comment: str


def process_publication(
    client: genai.Client, record: dict, article: dict, model: str
) -> Relevance:
    prompt = f"""\
You are a biomedical information extraction assistant.
You will be provided with a specific gene, a specific disease, and a \
scientific publication that may contain evidence for this gene-disease association.

Your task is to extract structured information only for the specified gene-disease pair \
, using only the information explicitly stated in the provided text.

For this association extract:
- mechanism: the explicitly stated mechanism of the disease. \
Only extract mechanisms that are explicitly and experimentally demonstrated in \
the text. Acceptable mechanisms must be supported by functional evidence or \
experimental assays described in the publication (e.g. "loss of function", "gain of function", "dominant negative"). \
Mechanisms must not be inferred from computational predictions (SIFT, PolyPhen), population data, \
phenotype correlations, or assumptions. Only use wording from the text or minimal paraphrase for clarity.
- mechanism_evidence: functional assay evidence directly supporting the \
stated mechanism of disease for the specified gene-disease pair. \
Only extract evidence if the publication explicitly describes a functional \
assay in which gene or protein function was directly measured and compared \
between normal and altered states (e.g. wild-type vs mutant) and which \
demonstrates the stated mechanism (e.g. loss of function, gain of function, dominant negative). \
Functional assay evidence can include experiments such as: \
protein expression or activity measurements; protein-protein or molecular interactions; \
cellular localization studies; rescue experiments, overexpression or \
knockdown/knockin studies; gene editing (CRISPR, morpholino) in cells or model organisms; \
assays in cell culture, mouse, zebrafish, Drosophila, or other systems; binding, \
enzymatic, or reporter assays \
Use the following keywords to identify functional assay evidence in the text: \
functional assays, mechanism, protein interaction, protein expression, interaction, expression, \
cells, cell culture, model organism, rescue, overexpression, gene editing, \
knockdown, knockin, crispr, morpholino, mouse, zebrafish, drosophila, \
cellular localisation, activity, binding.
- comment: a brief extractive justification that explicitly supports the extracted \
mechanism for the specified gene-disease association. \
The comment must only use wording present in the provided text or be a minimal paraphrase. \
Do not add new interpretations.

Only extract information if: \
- the gene and the disease are explicitly linked in the text \
- the mechanism or evidence is explicitly described

If the disease mentioned in the publication does not match the specified disease \
or if the gene-disease association is not explicitly supported \
return 'not found' for fields mechanism and mechanism_evidence.

Input:
Gene: {record['gene']}
Disease: {record['disease']}
Title: {article['title']}
Abstract: {article['abstract']}
Journal: {article['journal']}\
"""
    if article['fulltext']:
        prompt += "\nFull text: "+article['fulltext']

    response = client.models.generate_content(
        model=model,
        contents=prompt,
        config={
            "response_mime_type": "application/json",
            "response_schema": Relevance,
            "temperature": 0.2,
        },
    )

    return response.parsed


def get_article(pmid: int) -> dict | None:
    url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/article/MED/{pmid}?resultType=core&format=json"
    with urlopen(url) as res:
        data = json.loads(res.read())

    if data["hitCount"] == 0:
        return None

    obj = data["result"]

    full_text = None
    full_text_list = obj.get("fullTextIdList")
    if full_text_list:
        full_text_id = full_text_list.get("fullTextId")[0]

        url_full_text = f"https://www.ebi.ac.uk/europepmc/webservices/rest/{full_text_id}/fullTextXML"
        try:
            with urlopen(url_full_text) as response:
                xml_data = response.read()
                root = ET.fromstring(xml_data)
                title = root.findtext(".//article-title")
                abstract = " ".join(p.text.strip() for p in root.findall(".//abstract/p") if p.text)
                sections = []
                for sec in root.findall(".//body/*"):
                    sec_title = sec.findtext("title")
                    paragraphs = []

                    # Some section have the text directly there
                    for p in sec.findall("p"):
                        text = get_text_clean(p)
                        paragraphs.append(text.replace("\n", " "))
                    
                    # Other sections have sub-sections inside
                    for sub_section in sec.findall("sec"):
                        for p in sub_section.findall("p"):
                            sub_section_text = get_text_clean(p)
                            paragraphs.append(sub_section_text.replace("\n", " "))

                    sections.append((sec_title, paragraphs))

                full_text = ""
                for sec_title, paras in sections:
                    if sec_title:
                        full_text += "\n"+sec_title+"\n"
                    for p in paras:
                        full_text += " "+p
        except HTTPError as e:
            print(f"Error {e.code}: Full text not found for PMID {pmid}")

    title = obj.get("title")
    abstract = obj.get("abstractText")
    if not title or not abstract:
        return None

    parser = PlainTextExtractor()
    parser.feed(abstract)

    if "journalInfo" not in data["result"]:
        journal_info = "Unknown"
    else:
        journal_info = data["result"]["journalInfo"]["journal"]["title"]

    return {
        "title": title,
        "abstract": parser.get_text(),
        "journal": journal_info,
        "fulltext": full_text,
    }

def get_text_clean(element):
    """
    Extract all text inside an element.
    Keep text from <italic> and <bold>.
    Skip text inside <xref>, but preserve its tail text.
    """
    parts = []

    # Add leading text for the current element (if any)
    if element.text:
        parts.append(element.text)

    # Recursively process children
    for child in element:
        if child.tag == "xref":
            # Skip xref.text, but keep any tail text
            if child.tail:
                parts.append(child.tail)
        else:
            # Keep text and tail for all other tags (italic, bold, etc.)
            parts.append(get_text_clean(child))
            if child.tail:
                parts.append(child.tail)

    return "".join(parts)

class PlainTextExtractor(HTMLParser):
    def __init__(self):
        super().__init__()
        self.chunks = []

    def handle_data(self, data):
        if data.strip():
            self.chunks.append(data.strip())

    def get_text(self):
        return " ".join(self.chunks)


def download_g2p() -> list[dict]:
    """
    Method to download the G2P data from the API.
    It returns a list of all G2P records.
    """
    records = []

    # Download the DD data from the live API
    url = "https://www.ebi.ac.uk/gene2phenotype/api/panel/dd/download/"
    with urlopen(url) as res:
        data = res.read().decode("utf-8")

    f = StringIO(data)
    reader = csv.reader(f)
    it = iter(reader)
    keys = next(it)
    for values in it:
        obj = dict(zip(keys, values))

        publications = []
        for e in obj["additional mined publications"].split(";"):
            e = e.strip()
            if e:
                publications.append(
                    {
                        "id": int(e),
                        "title": None,
                        "abstract": None,
                        "journal": None,
                        "mechanism": None,
                        "mechanism_evidence": None,
                    }
                )
        for e in obj["publications"].split(";"):
            e = e.strip()
            if e:
                publications.append(
                    {
                        "id": int(e),
                        "title": None,
                        "abstract": None,
                        "journal": None,
                        "mechanism": None,
                        "mechanism_evidence": None,
                    }
                )

        record_to_append = {
                "id": obj["g2p id"],
                "gene": obj["gene symbol"],
                "disease": obj["disease name"],
                "confidence": obj["confidence"],
                "publications": publications,
            }

        records.append(record_to_append)

    return records


if __name__ == "__main__":
    main()
