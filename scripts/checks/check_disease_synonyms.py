#!/usr/bin/env python3

import os
import sys
import re
import argparse
import MySQLdb
import django
import configparser
from collections import defaultdict
from difflib import SequenceMatcher
from sentence_transformers import SentenceTransformer
from sklearn.metrics.pairwise import cosine_similarity

os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'gene2phenotype_app')
django.setup()
from gene2phenotype_app.utils import clean_string

"""
Description: Script to compare the similarity between the G2P diseases names and associated synonyms.
             Writes to the output file the diseases that have similarity score below the cutoff value.

Options
        --config : Config file with connection details to the G2P database (mandatory)
                File format is the following:
                        [g2p_database]
                        host = <>
                        port = <>
                        user = <>
                        password = <>
                        name = <>
        --cutoff : Score value to filter the results
        --output_file : Name of the output file

Requirements: Add the gene2phenotype_app project path to PYTHONPATH (export PYTHONPATH=_path_to_gene2phenotype_)
"""


def dump_g2p_diseases(db_host: str, db_port: int, db_name: str, user: str, password: str) -> dict[str, list]:
    """
    Queries the G2P database to dump all the diseases that have a disease synonym.

    Args:
        db_host (str): G2P database host name
        db_port (int): G2P database port number
        db_name (str): G2P database name
        user (str): user with read access
        password (str): password

    Returns:
        dict[str, list]: All disease names with associated synonyms
    """
    diseases = defaultdict(list)

    db = MySQLdb.connect(host=db_host, port=db_port, user=user, passwd=password, db=db_name)
    cursor = db.cursor()

    sql = """   SELECT d.name, s.synonym
                FROM disease d
                LEFT JOIN disease_synonym s ON s.disease_id = d.id
                WHERE s.synonym IS NOT NULL
          """

    cursor.execute(sql)
    data = cursor.fetchall()
    for row in data:
        disease_name = row[0]
        disease_synonym = row[1]
        diseases[disease_name].append(disease_synonym)

    return diseases


def compare_diseases(g2p_disease_synonyms : dict[str, list], cutoff: float, output_file: str) -> None:
    """
    Check how similar each G2P disease is with its associated synonyms.
    Write to a file the diseases with similarity score below the cutoff.

    Args:
        g2p_disease_synonyms (dict[str, list]): list of all diseases in G2P with associated synonyms list
        output_file (str): Output file name
        cutoff (float): return only diseases with similarity below this value
    """
    model = SentenceTransformer('pritamdeka/BioBERT-mnli-snli-scinli-scitail-mednli-stsb')

    with open(output_file, "w") as wr:
        wr.write(f"G2P disease name\tG2P disease synonym\tScore (SequenceMatcher)\tScore (BioBERT)\n")

        for disease_name, synonyms_list in g2p_disease_synonyms.items():
            # Clean the disease name
            new_disease_name = re.sub(".*\-related\s*", "", disease_name.lower()).strip()
            clean_disease_name = clean_string(new_disease_name)
            for synonym_name in synonyms_list:
                # Clean the synonym
                new_synonym_name = re.sub(".*\-related\s*", "", synonym_name.lower()).strip()
                clean_synonym_name = clean_string(new_synonym_name)

                # Calculate the string similarity
                score = SequenceMatcher(None, clean_disease_name, clean_synonym_name).ratio()

                embeddings = model.encode([clean_disease_name, clean_synonym_name])
                similarity = cosine_similarity([embeddings[0]], [embeddings[1]])

                if score < cutoff:
                    wr.write(f"{disease_name}\t{synonym_name}\t{score}\t{similarity[0][0]}\n")


def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--config", required=True, help="Config file with details to the G2P database")
    parser.add_argument("--cutoff", default=0.5, help="Return scores below this cutoff")
    parser.add_argument("--output_file", default="report_disease_scores.txt", help="Output file name")

    args = parser.parse_args()
    config_file = args.config
    cutoff = args.cutoff
    output_file = args.output_file

    # Load the config file
    config = configparser.ConfigParser()
    config.read(config_file)

    try:
        g2p_config = config['g2p_database']
    except KeyError:
        sys.exit("ERROR: 'g2p_database' missing from config file")
    else:
        db_host = g2p_config['host']
        db_port = g2p_config['port']
        db_name = g2p_config['name']
        user = g2p_config['user']
        password = g2p_config['password']

    g2p_disease_synonyms = dump_g2p_diseases(db_host, int(db_port), db_name, user, password)
    compare_diseases(g2p_disease_synonyms, float(cutoff), output_file)

if __name__ == '__main__':
    main()