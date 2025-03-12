#!/usr/bin/env python3

import os.path
import sys
import argparse
import MySQLdb
import requests
import configparser

lgd_disease_to_update = []

def dump_data(db_host, db_port, db_name, user, password):
    gene_records = {} # key = gene symbol; value = list of records
    diseases = {} # key = disease name; value = list of record ids

    sql = """   SELECT l.name, g2p.stable_id, d.name, d.id, a.value, a1.value, lgd.id
                FROM locus_genotype_disease lgd
                LEFT JOIN locus l ON l.id = lgd.locus_id
                LEFT JOIN g2p_stableid g2p ON g2p.id = lgd.stable_id
                LEFT JOIN disease d ON d.id = lgd.disease_id
                LEFT JOIN attrib a ON a.id = lgd.genotype_id
                LEFT JOIN attrib a1 ON a1.id = lgd.confidence_id
                ORDER BY l.name
          """
    
    sql_disease = """
                        SELECT d.name, d.id, lgd.id
                        FROM disease d
                        LEFT JOIN locus_genotype_disease lgd ON lgd.disease_id = d.id
                        ORDER BY d.name
                  """

    db = MySQLdb.connect(host=db_host, port=db_port, user=user, passwd=password, db=db_name)
    cursor = db.cursor()
    cursor.execute(sql)
    data = cursor.fetchall()
    for row in data:
        if row[0] not in gene_records:
            gene_records[row[0]] = [{
                "stable_id": row[1],
                "disease_name": row[2],
                "disease_id": row[3],
                "genotype": row[4],
                "confidence": row[5],
                "record_id": row[6]
            }]
        else:
            gene_records[row[0]].append({
                "stable_id": row[1],
                "disease_name": row[2],
                "disease_id": row[3],
                "genotype": row[4],
                "confidence": row[5],
                "record_id": row[6]
            })

    cursor.execute(sql_disease)
    data_disease = cursor.fetchall()
    for row in data_disease:
        if row[0] not in diseases:
            diseases[row[0]] = {
                "disease_id": row[1],
                "records": [row[2]]
            }
        else:
            diseases[row[0]]["records"].append(row[2])

    db.close()

    return gene_records, diseases

def read_file(file, gene_records, diseases, dryrun):
    file_output = "gene_disease_found_in_g2p.txt"
    file_not_updated = "diseases_not_updated.txt"
    diseases_to_update = []

    with open(file, "r") as fh, open(file_output, "w") as wr, open(file_not_updated, "w") as wr_diseases:
        wr.write("gene symbol\tdisease name\tdisease name formatted\tdiseases found in G2P linked to gene\n")

        # Header:
        # gene symbol, disease name, disease name formatted, Updated
        # Other columns are ignored
        for line in fh:
            if not line.startswith("gene symbol"):
                data = line.split("\t")
                gene_symbol = data[0]
                current_disease = data[1] # before searching for the disease we have to check if 'gene-related' is used
                new_disease = data[2]
                is_updated = data[3]

                if is_updated == "Yes":
                    continue

                # Check if the record can be found in G2P
                try:
                    db_data = gene_records[gene_symbol]
                except KeyError:
                    print(f"WARNING: {gene_symbol} not found in G2P")
                else:
                    to_update = 0
                    list_disease = []

                    if dryrun:
                        print(f"\n{gene_symbol}; {current_disease} (New disease: {new_disease})")

                    for record in db_data:
                        # In the old system most of the disease names don't include the '<gene>-related' in the name
                        # to compare disease names we have to compare the end of the name
                        if ((record["disease_name"].lower() == current_disease.lower() or record["disease_name"].endswith("related "+current_disease)
                             or record["disease_name"].endswith("associated "+current_disease))
                            and record["disease_name"] != new_disease):
                            to_update = 1
                            # Check if new disease name is already in the db
                            if new_disease in diseases:
                                # If the new disease name already exists then add it to replace existing name
                                # Before running the updating, the API checks if records are not the same
                                print(f"Update LGD records: replace disease_id {record['disease_id']} by {diseases[new_disease]['disease_id']}")
                                lgd_disease_to_update.append(
                                    {
                                        "disease_id": record["disease_id"],
                                        "new_disease_id": diseases[new_disease]["disease_id"]
                                    }
                                )

                                if dryrun:
                                    print(f"Update disease in lgd table -> replace disease_id: {record['disease_id']} by {diseases[new_disease]['disease_id']}")

                            elif gene_symbol not in new_disease:
                                # This should not happen
                                # Action: print to file 'diseases_not_updated.txt'
                                wr_diseases.write("Gene not found in new disease name\t"+line)
                            else:
                                # Create list of diseases to update
                                diseases_to_update.append(
                                    {
                                        "id": record["disease_id"],
                                        "name": new_disease
                                    }
                                )
                                if dryrun:
                                    print(f"Update disease name -> disease_id: {record['disease_id']}; current name: {current_disease}; new name: {new_disease}")

                        else:
                            list_disease.append(record["disease_name"])

                    # Print to file genes that won't have disease updates
                    if not to_update:
                        wr.write(f"{gene_symbol}\t{current_disease}\t{new_disease}\t{list_disease}\n")

    return diseases_to_update

def update_diseases(diseases_to_update, api_username, api_password, api_url):
    """
        Method to update the disease names.
        This method calls the gene2phenotype api and updates the database defined in gene2phenotype_project/config.ini
    """
    disease_url = "update/diseases/"
    lgd_disease_url = "lgd_disease_updates/"
    login_url = "login/"

    data = {
        "username": api_username,
        "password": api_password
    }

    response = requests.post(api_url + login_url, json=data)
    if response.status_code == 200:
        try:
            response_update = requests.post(api_url + disease_url, json=diseases_to_update, cookies=response.cookies)
            if response_update.status_code == 200:
                response_json = response_update.json()
                print("Diseases updated successfully:", response_json)
                # Some diseases were probably not updated: if they are already in the db
                # For those we can update the disease_id in lgd to point to the existing disease
                if "errors" in response_json:
                    for error in response_json["errors"]:
                        print(f"Update LGD records: replace disease_id {error['id']} by {error['existing_id']}")
                        lgd_disease_to_update.append(
                            {
                                "disease_id": error['id'],
                                "new_disease_id": error['existing_id']
                            }
                        )
                    response_update_lgd = requests.post(api_url + lgd_disease_url, json=lgd_disease_to_update, cookies=response.cookies)
                    if response_update_lgd.status_code == 200:
                        print("LGD records updated successfully:", response_update_lgd.json())
                    else:
                        print("Failed to update LGD records:", response_update_lgd.status_code, response_update_lgd.json())
            else:
                print("Failed to update diseases:", response_update.status_code, response_update.json())
        except Exception as e:
            print("Error:", e)
    else:
        print("Error: cannot login into G2P")


def main():
    """
        Params:
                --config : Config file name containing the database and API connection info (mandatory)
                        File format is the following: 
                            [database]
                            host = <>
                            port = <>
                            user = <>
                            password = <>
                            name = <>

                            [api]
                            api_url = <>

                --file : Tab delimited file with all diseases to be updated (mandatory)
                    File format is the following:
                        gene symbol\tdisease name\tdisease name formatted\tUpdated
                
                --api_username: Username to connect to the G2P API (mandatory)
                --api_password: Password to connect to the G2P API (mandatory)
                --dryrun: Test script without running the updates (default: 0)
    """

    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--config", required=True, help="Config file")
    parser.add_argument("--file", required=True, help="Tab delimited file with all diseases to be updated")
    parser.add_argument("--api_username", required=True, help="Username to connect to the G2P API")
    parser.add_argument("--api_password", required=True, help="Password to connect to the G2P API")
    parser.add_argument("--dryrun", required=False, default=0, help="Option to test update")
    args = parser.parse_args()

    file = args.file
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

    print("Dump data from G2P...")
    gene_records, diseases = dump_data(db_host, int(db_port), db_name, user, password)
    print("Dump data from G2P... done\n")

    if os.path.isfile(file):
        expected_columns = ['gene symbol', 'disease name', 'disease name formatted', 'Updated']

        with open(file, "r", encoding="utf-8") as fh:
            header = fh.readline().strip().split("\t")
            if header[:4] != expected_columns:
                sys.exit(f"Error: File format is incorrect. Found: {header[:4]}; Expected: {expected_columns}")

        print("Parsing diseases to update...")
        diseases_to_update = read_file(file, gene_records, diseases, dryrun)
        print("Parsing diseases to update... done\n")

        if not dryrun:
            print("Updating disease names...")
            update_diseases(diseases_to_update, api_username, api_password, api_url)
            print("Updating disease names... done\n")

    else:
        print(f"Input file is invalid '{file}'")


if __name__ == '__main__':
    main()