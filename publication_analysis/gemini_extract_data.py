#!/usr/bin/env python3

import argparse
import configparser
import csv
import json
import sys
import time
from html.parser import HTMLParser
from pathlib import Path
from urllib.request import urlopen
from urllib.error import HTTPError
import xml.etree.ElementTree as ET

from google import genai
from google.genai.types import HttpOptions
from google.oauth2 import service_account
from pydantic import BaseModel


"""
Script extract data from a publication (PMID) using Google Gemini.
Input files may be CSV (.csv) or TSV, with PMIDs in the first column.

Usage:
    $ python gemini_extract_data.py --config config.ini --limit 300 --output g2p_records_mined.json g2p_records.tsv

The config.ini file should contain the following entries:
    [project_config]
    key_file = path/to/key_file.json
    project = your-project-id
    location = europe-west2
    model = gemini-2.5-flash
"""


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("infile", type=Path, help="Input CSV (.csv) or TSV file with PMIDs (first column)")
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Output pretty-printed JSON file"
    )
    parser.add_argument(
        "--config",
        type=Path,
        required=True,
        help="Config file with Google settings(key file, project name, model, etc.)",
    )
    parser.add_argument(
        "-l",
        "--limit",
        type=int,
        default=0,
        metavar="N",
        help="Process N publications and exit (default: all)",
    )
    parser.set_defaults(func=run_process)

    args = parser.parse_args()
    args.func(args)


def run_process(args):
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
        location=config["location"],
        credentials=credentials,
        http_options=HttpOptions(api_version="v1")
    )

    existing_pmids = load_existing_pmids(args.output)
    if existing_pmids:
        print(
            f"Skipping {len(existing_pmids)} PMIDs already present in {args.output}",
            file=sys.stderr,
        )

    pmids_to_analyse = {}
    delimiter = "," if args.infile.suffix.lower() == ".csv" else "\t"
    with args.infile.open("rt", encoding="utf-8-sig", newline="") as fh:
        for fields in csv.reader(fh, delimiter=delimiter):
            if not fields:
                continue
            pmid = fields[0].strip()
            if not pmid or pmid.upper() == "PMID" or pmid in existing_pmids:
                continue
            if pmid not in pmids_to_analyse:
                pmids_to_analyse[pmid] = {
                    "pmid": pmid,
                    "results": [],
                    "status": None,
                }

    try:
        done = 0

        for pub in pmids_to_analyse.values():
            article = get_article(pub["pmid"])
            if article is None:
                pub["status"] = "incomplete"
                append_json_record(args.output, pub)
                existing_pmids.add(pub["pmid"])
                done += 1
                if args.limit > 0 and done == args.limit:
                    return
                continue

            pub.update(**article)
            pub["content_source"] = "fulltext" if pub.get("fulltext") else "abstract"

            output = process_publication(client, pub, config["model"])
            pub["results"] = [result.model_dump() for result in output.results]
            pub["status"] = "complete"

            if done:
                print("", file=sys.stderr)

            print(
                f"PMID        : {pub['pmid']}",
                file=sys.stderr,
            )
            print(f"Article  : {pub['title']}", file=sys.stderr)
            print(f"Source   : {pub['content_source']}", file=sys.stderr)
            print(f"Results  : {len(pub['results'])}", file=sys.stderr)

            append_json_record(args.output, pub)
            existing_pmids.add(pub["pmid"])

            done += 1
            if args.limit > 0 and done == args.limit:
                return

    finally:
        print(f"Processed {done} publications", file=sys.stderr)


def load_json_key(key_file):
    credentials = service_account.Credentials.from_service_account_file(key_file).with_scopes(["https://www.googleapis.com/auth/cloud-platform"])
    return credentials


def load_existing_pmids(output_path: Path) -> set[str]:
    if not output_path.exists() or output_path.stat().st_size == 0:
        return set()

    try:
        with output_path.open("rt", encoding="utf-8") as fh:
            data = json.load(fh)
    except json.JSONDecodeError as exc:
        raise SystemExit(
            f"Existing output file is not valid JSON: {output_path}. "
            "Fix or remove it before resuming."
        ) from exc

    if not isinstance(data, list):
        raise SystemExit(f"Existing output file must contain a JSON list: {output_path}")

    return {str(item.get("pmid")) for item in data if isinstance(item, dict) and item.get("pmid")}


def append_json_record(output_path: Path, record: dict) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_record = {key: value for key, value in record.items() if key != "fulltext"}

    if not output_path.exists() or output_path.stat().st_size == 0:
        with output_path.open("wt", encoding="utf-8") as fh:
            fh.write("[\n")
            json.dump(output_record, fh, indent=2, ensure_ascii=False)
            fh.write("\n]\n")
        return

    with output_path.open("rb+") as fh:
        fh.seek(0, 2)
        pos = fh.tell() - 1

        while pos >= 0:
            fh.seek(pos)
            char = fh.read(1)
            if char not in b" \t\r\n":
                break
            pos -= 1

        if pos < 0:
            raise SystemExit(f"Existing output file is malformed: {output_path}")

        fh.seek(pos)
        if fh.read(1) != b"]":
            raise SystemExit(
                f"Existing output file must end with a closing JSON list bracket: {output_path}"
            )

        previous_pos = pos - 1
        while previous_pos >= 0:
            fh.seek(previous_pos)
            char = fh.read(1)
            if char not in b" \t\r\n":
                break
            previous_pos -= 1

        if previous_pos < 0:
            raise SystemExit(f"Existing output file is malformed: {output_path}")

        has_records = char != b"["
        fh.seek(pos)
        fh.truncate()
        if has_records:
            fh.write(b",\n")
        fh.write(json.dumps(output_record, indent=2, ensure_ascii=False).encode("utf-8"))
        fh.write(b"\n]\n")


class GeneDiseaseMechanism(BaseModel):
    gene: str
    disease: str
    mode_of_inheritance: str
    mechanism: str
    mechanism_evidence: list[str]
    number_of_families: int | None
    number_of_affected_individuals: int | None
    comment: str


class RelevanceResponse(BaseModel):
    results: list[GeneDiseaseMechanism]


def process_publication(
    client: genai.Client, article: dict, model: str
) -> RelevanceResponse:
    prompt = f"""\
You are a biomedical information extraction assistant.
You will be provided with a scientific publication that may contain evidence \
for a gene-disease association.

Your task is to return a list of structured gene-disease mechanism records, \
using only the information explicitly stated in the provided text.

For each supported gene-disease association extract:
- gene: the gene symbol (e.g. "BRCA1").
- disease: the associated disease name (e.g. "breast cancer").
- mode_of_inheritance: the explicitly stated mode of inheritance or allelic requirement for the \
gene-disease association (e.g. "autosomal dominant", "autosomal recessive", "biallelic", \
"monoallelic"). Only extract what is explicitly stated in the provided text. If mode of \
inheritance is not available, set this field to an empty string.
- mechanism: the explicitly stated mechanism of the disease. \
Only extract mechanisms that are explicitly and experimentally demonstrated in \
the text. Acceptable mechanisms must be supported by functional evidence or \
experimental assays described in the publication (e.g. "loss of function", "gain of function", "dominant negative"). \
Mechanisms must not be inferred from computational predictions (SIFT, PolyPhen), population data, \
phenotype correlations or assumptions. Only use wording from the text or minimal paraphrase for clarity.
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
enzymatic, or reporter assays. \
Use the following keywords to identify functional assay evidence in the text: \
functional assays, mechanism, protein interaction, protein expression, interaction, expression, \
cells, cell culture, model organism, rescue, overexpression, gene editing, \
knockdown, knockin, crispr, morpholino, mouse, zebrafish, drosophila, \
cellular localisation, activity, binding.
- number_of_families: the explicit number of families reported for the specified \
gene-disease association in the current publication. Only extract a numeric count \
if it is directly stated in the provided text and clearly tied to this specific \
gene-disease association. Do not infer from tables, cohorts, pedigrees, or prior \
studies unless the count is explicitly stated for this association in the current \
article. If unavailable, vague (e.g. "several families"), or only stated for the \
paper more generally rather than this association, set to null.
- number_of_affected_individuals: the explicit number of affected individuals \
reported for the specified gene-disease association in the current publication. \
Only extract a numeric count if it is directly stated in the provided text and \
clearly tied to this specific gene-disease association. Do not infer by counting \
cases, pedigrees, or variants, and do not use counts reported only in previous \
studies. If unavailable, vague (e.g. "multiple affected individuals"), or only \
stated for the paper more generally rather than this association, set to null.
- comment: a brief extractive justification that explicitly supports the extracted \
mechanism for the specified gene-disease association. \
The comment must only use wording present in the provided text or be a minimal paraphrase. \
Do not add new interpretations.

Return JSON with a top-level field `results`, where `results` is a list of objects. \
Each object must contain exactly these fields: `gene`, `disease`, \
`mode_of_inheritance`, `mechanism`, `mechanism_evidence`, \
`number_of_families`, `number_of_affected_individuals`, and `comment`.

If the publication supports a gene-disease association but does not provide valid \
mechanism evidence, still return the gene-disease pair in `results`. In that case:
- set `mode_of_inheritance` to an empty string if it is not available
- set `mechanism` to an empty string
- set `mechanism_evidence` to an empty list
- set `number_of_families` to null if it is not available
- set `number_of_affected_individuals` to null if it is not available
- set `comment` to an empty string
Do not omit a supported gene-disease pair solely because mechanism information \
is unavailable.

Do not extract mechanism or mechanism_evidence if: \
- the evidence is based on functional studies from previous studies or publications \
- the evidence is based on computational predictions, population data, phenotype correlations or assumptions

Treat review-style, background, and literature-summary statements as non-primary evidence. \
Do not use them as mechanism evidence unless the current paper itself reports the experiment. \
In particular, do not treat phrases such as "has been shown", "previous studies", \
"previously reported", "it was reported", "in earlier work", "others showed", \
"is known to", "reviewed in", or citation-supported summary statements as primary evidence \
for mechanism or mechanism_evidence unless the current article explicitly describes the assay \
performed in this study.

Input:
Title: {article['title']}
Abstract: {article['abstract']}
"""
    if article['fulltext']:
        prompt += "\nFull text: "+article['fulltext']

    response = client.models.generate_content(
        model=model,
        contents=prompt,
        config={
            "response_mime_type": "application/json",
            "response_schema": RelevanceResponse,
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


if __name__ == "__main__":
    main()
