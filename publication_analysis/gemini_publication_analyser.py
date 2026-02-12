#!/usr/bin/env python3

import argparse
import csv
import configparser
import json
import sys
import time
from enum import Enum
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
Script to download G2P records and associated mined publications and to assess the
relevance of each publication to the G2P record using Google Vertex AI Gemini models.

Usage:
    1) Download G2P records and associated mined publications
    $ python gemini_publication_analyser.py init g2p_records.json

    2) Assess the relevance for 300 publications
    $ python gemini_publication_analyser.py process --config config.ini --limit 300 g2p_records.json

The config.ini file should contain the following entries:
    [project_config]
    key_file = path/to/key_file.json
    project = your-project-id
    location = europe-west2
    model = gemini-2.5-flash
"""


def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(required=True)

    p1 = subparsers.add_parser("init")
    p1.add_argument("output", type=Path, help="Output JSON file")
    p1.set_defaults(func=run_download)

    p2 = subparsers.add_parser("process")
    p2.add_argument("infile", type=Path, help="Input JSON file")
    p2.add_argument(
        "--config",
        type=Path,
        required=True,
        help="Config file with Google Vertex AI settings(key file, project name, model, etc.)",
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
        "--confidence",
        type=str,
        default=None,
        help="Process only records with specific confidence value (default: all)",
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

    try:
        done = 0
        
        if args.rpm > 0:
            secs = 60.0 / args.rpm
        else:
            secs = 0

        for record in records:
            if args.confidence:
                if record["confidence"] != args.confidence:
                    continue

            for pub in record["publications"]:
                if pub["status"] is not None:
                    continue

                article = get_article(pub["id"])
                if article is None:
                    pub["status"] = "incomplete"
                    continue

                pub.update(**article)

                relevance = process_publication(client, record, pub, config["model"])
                pub["status"] = relevance.label.value
                pub["comment"] = relevance.comment
                pub["ai_model"] = config["model"]

                if done:
                    print("", file=sys.stderr)

                if "mechanism" in record:
                    mechanism_value = record['mechanism']
                else:
                    mechanism_value = "undetermined"

                print(
                    f"G2P      : {record['id']}, {record['gene']}, {mechanism_value}, {record['disease']}, {record['confidence']}",
                    file=sys.stderr,
                )
                print(f"Article  : {pub['title']}", file=sys.stderr)
                print(f"Relevance: {pub['status']}", file=sys.stderr)
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


class Label(Enum):
    HIGH = "high"
    MEDIUM = "medium"
    LOW = "low"
    DISPUTED = "disputed"


class Relevance(BaseModel):
    label: Label
    comment: str


def process_publication(
    client: genai.Client, record: dict, article: dict, model: str
) -> Relevance:
    prompt = f"""\
You are assessing whether a scientific publication is relevant \
to a Gene2Phenotype record defined by a gene and a disease. 
Relevance means the publication provides or discusses evidence \
that this gene is causally or mechanistically linked to this \
disease in humans or relevant models (e.g. mammalian or functional \
models recapitulating the human phenotype).

Output one of four labels:
- high: the article directly supports or reports an association \
between the specified gene and the specified disease (same disease, \
not just related systems).
- medium: the article discusses the specified gene or disease in a \
mechanistically relevant or closely related context, but without \
demonstrating a direct association between the two. \
The disease context should still be similar \
(e.g. same system or phenotype family).
- low: the article discusses the gene or disease in an unrelated context, \
or links the gene to a different disease than the one specified or focuses on \
a different gene.
- disputed: The article provides evidence that contradicts or disproves \
an association between the specified gene and the specified disease.

Then provide one short reason.

NEVER assign "high" relevance unless \
the publication provides evidence directly linking the gene to the \
specific disease named in the record.
If the article discusses a different disease caused by the same gene: \
- If the diseases share overlapping molecular mechanisms and phenotypes, assign "medium"; \
- If they do not, assign "low".
Consider whether the molecular mechanism described in the publication \
(e.g. gain or loss of function) matches the mechanism in the record (if available) \
when assessing similarity.
If the publication or the record do not mention a molecular mechanism, base your \
decision on the gene–disease association itself (e.g. clinical or genetic evidence). \
Do not lower relevance solely because the mechanism is unspecified.
If the publication discusses multiple genes or structural variants involving \
the specified gene then assign "low" relevance.

Input:
Gene: {record['gene']}
Previous gene symbols: {record['previous_gene_symbols']}
Disease: {record['disease']}
Title: {article['title']}
Abstract: {article['abstract']}
Journal: {article['journal']}\
"""
    if "mechanism" in record:
        prompt += "\nMolecular mechanism: "+record['mechanism']

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
                        "status": None,
                        "comment": None,
                    }
                )

        record_to_append = {
                "id": obj["g2p id"],
                "gene": obj["gene symbol"],
                "previous_gene_symbols": obj["previous gene symbols"],
                "disease": obj["disease name"],
                "confidence": obj["confidence"],
                "publications": publications,
            }
        if obj["molecular mechanism"] != "undetermined":
            record_to_append["mechanism"] = obj["molecular mechanism"]

        records.append(record_to_append)

    return records


if __name__ == "__main__":
    main()
