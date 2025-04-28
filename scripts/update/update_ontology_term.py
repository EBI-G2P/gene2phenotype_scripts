#!/usr/bin/env python3

import sys
import argparse
import MySQLdb
import requests
import configparser
import re


def dump_ontology_data(db_host: str, db_port: int, db_name: str, user: str, password: str) -> dict[str, dict]:
    """
    Method to fetch all the OMIM and Mondo ontology terms stored in G2P.

    Args:
        db_host (str): G2P database host name
        db_port (int): G2P database port number
        db_name (str): G2P database name
        user (str): user with read access
        password (str): password

    Returns:
        dict[str, dict]: dictionary with ontology accession and its corresponding data
    """

    ontology_records = {}

    sql = """   SELECT o.accession, o.term, o.description
                FROM ontology_term o
                LEFT JOIN source s ON s.id = o.source_id
                WHERE s.name = 'Mondo' OR s.name = 'OMIM'
          """

    db = MySQLdb.connect(host=db_host, port=db_port, user=user, passwd=password, db=db_name)
    cursor = db.cursor()
    cursor.execute(sql)
    data = cursor.fetchall()
    for row in data:
        try:
            ontology_records[row[0]]
        except KeyError:
            ontology_records[row[0]] = {
                "term": row[1],
                "description": row[2]
            }
        else:
            print(f"WARNING: duplicated ontology term '{row[0]}'")

    db.close()

    return ontology_records


def analyse_terms(ontology_records: dict[str, dict]) -> tuple[dict[str, dict], list[str]]:
    """
    Method to analyse which ontology terms have to be updated or removed.
    Records to update are:
        - ontology terms with accession equal to the term

    Args:
        ontology_records (dict[str, dict]): dictionary with all ontology terms

    Returns:
        dict[str, dict]: dictionary with ontology terms to be updated
        list[str]: list of ontology IDs that are obsolete
    """    
    records_to_update = {}
    records_to_delete = []

    for accession in ontology_records:
        # Checks if the term name is the same as the accession
        # For this scenario we have to query the Mondo API to get the
        # correct term
        if accession == ontology_records[accession]["term"]:
            if accession.startswith("MONDO"):
                new_term, new_description = get_mondo(accession)
            else:
                new_term, new_description = get_omim(accession)
            
            # If the ontology ID has a term then add it to the list of records to be updated
            if new_term:
                records_to_update[accession] = {
                    "term": new_term,
                    "description": new_description
                }
            else:
                print(f"WARNING: invalid ontology ID {accession}")
                records_to_delete.append(accession)

    return records_to_update, records_to_delete


def get_mondo(id: str) -> tuple[str, str]:
    """
    Fetch the disease term and description for a specific Mondo ID.

    Args:
        id (str): Mondo ID

    Returns:
        tuple[str, str]: term and description
    """

    url = f"https://www.ebi.ac.uk/ols4/api/search?q={id}&ontology=mondo&exact=1"

    term = None
    description = None

    r = requests.get(url, headers={ "Content-Type" : "application/json"})

    if not r.ok:
        return None

    decoded = r.json()
    if decoded["response"]["numFound"] != 0:
        if len(decoded['response']['docs'][0]['description']) > 0:
            description = decoded['response']['docs'][0]['description'][0]
        term = decoded['response']['docs'][0]['label']

    return term, description


def get_omim(id: str) -> tuple[str, str]:
    """
    Get OMIM data from OLS API.

    Args:
        id (str): OMIM ID

    Returns:
        tuple[str, str]: ontology term and description
    """    
    disease = None
    description = None

    url = f"https://www.ebi.ac.uk/ols4/api/search?q={id}&ontology=cco"

    r = requests.get(url, headers={ "Content-Type" : "application/json"})

    if not r.ok:
        return disease, description

    decoded = r.json()

    if decoded["response"]["numFound"] != 0:
        source_name = decoded["response"]["docs"][0]["iri"]
        # OLS API can return data from different sources
        # Check only for OMIM data
        if source_name.find("/omim/") != -1:
            disease = decoded['response']['docs'][0]['label']

            if len(decoded['response']['docs'][0]['description']) > 0:
                description = decoded['response']['docs'][0]['description'][0]

    return disease, description


def update_ontologies(records_to_update: dict[str, dict], api_username: str, api_password: str, api_url: str) -> None:
    """
    Method that calls the G2P API to update all the ontology terms.

    Args:
        records_to_update (dict[str, dict]): dictionary with ontology terms to be updated
        api_username (str): G2P API username
        api_password (str): G2P API password
        api_url (str): G2P API URL
    """

    ontology_url = "update/ontology_terms/"
    login_url = "login/"

    data = {
        "username": api_username,
        "password": api_password
    }

    response = requests.post(api_url + login_url, json=data)
    if response.status_code == 200:
        try:
            response_update = requests.post(api_url + ontology_url, json=records_to_update, cookies=response.cookies)
            if response_update.status_code == 200:
                response_json = response_update.json()
                print("Ontologies updated successfully:", response_json)
                if "errors" in response_json:
                    for error in response_json["errors"]:
                        print(f"ERROR: {error}")
            else:
                print("Failed to update ontologies:", response_update.status_code, response_update.json())
        except Exception as e:
            print("Error:", e)
    else:
        print("Error: cannot login into G2P")


def delete_ontologies(records_to_delete: list[str], api_username: str, api_password: str, api_url: str) -> None:
    """
    Method that calls the G2P API to delete a list of ontology accession.

    Args:
        records_to_delete (list[str]): list of ontology accessions to delete
        api_username (str): G2P API username
        api_password (str): G2P API password
        api_url (str): G2P API URL
    """

    ontology_url = "update/ontology_terms/"
    login_url = "login/"

    data = {
        "username": api_username,
        "password": api_password
    }

    # Login
    response = requests.post(api_url + login_url, json=data)
    if response.status_code == 200:
        try:
            # Delete ontology terms
            response_delete = requests.delete(api_url + ontology_url, json=records_to_delete, cookies=response.cookies)
            if response_delete.status_code == 200:
                response_json = response_delete.json()
                print("Ontologies delete successfully:", response_json)
                if "errors" in response_json:
                    for error in response_json["errors"]:
                        print(f"ERROR: {error}")
            else:
                print("Failed to delete ontologies:", response_delete.status_code, response_delete.json())
        except Exception as e:
            print("Error:", e)
    else:
        print("Error: cannot login into G2P")


def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--config", required=True, help="Config file")
    parser.add_argument("--api_username", required=True, help="Username to connect to the G2P API")
    parser.add_argument("--api_password", required=True, help="Password to connect to the G2P API")
    parser.add_argument("--dryrun", required=False, default=0, help="Option to test update")
    args = parser.parse_args()

    config_file = args.config
    dryrun = args.dryrun
    api_username = args.api_username
    api_password = args.api_password

    # Load the config file
    config = configparser.ConfigParser()
    config.read(config_file)

    db_host = config['database']['host']
    db_port = config['database']['port']
    db_name = config['database']['name']
    user = config['database']['user']
    password = config['database']['password']
    api_url = config['api']['api_url']

    print("Dump ontology data from G2P...")
    ontology_records = dump_ontology_data(db_host, int(db_port), db_name, user, password)
    print("Dump ontology data from G2P... done\n")

    records_to_update, records_to_delete = analyse_terms(ontology_records)

    if not dryrun:
        print("Updating ontologies...")
        update_ontologies(records_to_update, api_username, api_password, api_url)
        print("Updating ontologies... done\n")
        print("Deleting ontologies...")
        delete_ontologies(records_to_delete, api_username, api_password, api_url)
        print("Deleting ontologies... done\n")
    else:
        print("Updates to run:")
        for record in records_to_update:
            print(f"{record}: {records_to_update[record]['term']}")
        print("\nOntologies to delete:", records_to_delete)

if __name__ == '__main__':
    main()