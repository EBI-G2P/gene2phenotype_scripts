#!/usr/bin/env python3

import os.path
import sys
import argparse
import MySQLdb
import requests
import configparser

"""
    Script to update disease names.
    If the new disease name already exists in the db, then it updates
    the disease id in locus_genotype_disease with the existing disease id.
    It prints an info message if the disease is linked to more than one gene.

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
                    g2p id\tgene symbol\tdisease name\tdisease name formatted\tUpdated

            --api_username: Username to connect to the G2P API (mandatory)
            --api_password: Password to connect to the G2P API (mandatory)
            --dryrun: Test script without running the updates (default: 0)
"""

# List of disease ids to be replaced by an existing disease id
# The update is done in the locus_genotype_disease (replace disease id)
lgd_disease_to_update = []
# Save a list of unique diseases from the input file - to help keeping track of the updates
unique_diseases_from_input = set()


def dump_data(
    db_host: str, db_port: int, db_name: str, user: str, password: str
) -> tuple[dict[str, list], dict[str, dict]]:
    """
    Dumps the current data from the G2P database.

    Args:
        db_host (str): hostname
        db_port (int): port
        db_name (str): G2P database name
        user (str): username
        password (str): password

    Returns:
        tuple[dict[str, list], dict[str, dict]]: dict of genes and associated diseases,
        and a dict of diseases by name and associated records
    """
    gene_records = {}  # key = gene symbol; value = list of records
    diseases = {}  # key = disease name; value = list of record ids

    sql = """   SELECT l.name, g2p.stable_id, d.name, d.id, a.value, a1.value, lgd.id
                FROM locus_genotype_disease lgd
                LEFT JOIN locus l ON l.id = lgd.locus_id
                LEFT JOIN g2p_stableid g2p ON g2p.id = lgd.stable_id
                LEFT JOIN disease d ON d.id = lgd.disease_id
                LEFT JOIN attrib a ON a.id = lgd.genotype_id
                LEFT JOIN attrib a1 ON a1.id = lgd.confidence_id
                ORDER BY l.name
          """

    sql_disease = """ SELECT d.name, d.id, lgd.id, g.stable_id
                      FROM disease d
                      LEFT JOIN locus_genotype_disease lgd ON lgd.disease_id = d.id
                      LEFT JOIN g2p_stableid g ON g.id = lgd.stable_id
                      ORDER BY d.name
                  """

    db = MySQLdb.connect(
        host=db_host, port=db_port, user=user, passwd=password, db=db_name
    )
    cursor = db.cursor()
    cursor.execute(sql)
    data = cursor.fetchall()
    for row in data:
        if row[0] not in gene_records:
            gene_records[row[0]] = [
                {
                    "stable_id": row[1],
                    "disease_name": row[2],
                    "disease_id": row[3],
                    "genotype": row[4],
                    "confidence": row[5],
                    "record_id": row[6],
                }
            ]
        else:
            gene_records[row[0]].append(
                {
                    "stable_id": row[1],
                    "disease_name": row[2],
                    "disease_id": row[3],
                    "genotype": row[4],
                    "confidence": row[5],
                    "record_id": row[6],
                }
            )

    cursor.execute(sql_disease)
    data_disease = cursor.fetchall()
    for row in data_disease:
        if row[0] not in diseases:
            diseases[row[0]] = {
                "disease_id": row[1],
                "records": [{"lgd_id": row[2], "stable_id": row[3]}],
            }
        else:
            diseases[row[0]]["records"].append({"lgd_id": row[2], "stable_id": row[3]})

    db.close()

    return gene_records, diseases


def read_file(
    file: str, gene_records: dict[str, list], diseases: dict[str, dict], dryrun: int
) -> list:
    """
    Reads the diseases from the input file, for each disease checks if it can be udpdated.

    Args:
        file (str): tab delimited input file
        gene_records (dict[str, list]): dictionary of genes and associated list of diseases
        diseases (dict[str, dict]): dictionary of diseases by name and associated records
        dryrun (int): run script in update mode (default: 0)

    Returns:
        list: list of diseases to be updated
    """
    file_output = "gene_disease_found_in_g2p.txt"
    file_not_updated = "diseases_not_updated.txt"
    diseases_to_update = []

    with (
        open(file, "r") as fh,
        open(file_output, "w") as wr,
        open(file_not_updated, "w") as wr_diseases,
    ):
        wr.write(
            "g2p id\tgene symbol\tdisease name\tdisease name formatted\tdiseases found in G2P linked to gene\n"
        )

        # Header:
        # g2p id, gene symbol, disease name, disease name formatted, Updated
        # Other columns are ignored
        for line in fh:
            if not line.startswith("g2p id"):
                data = line.strip().split("\t")
                g2p_id = data[0].strip()
                gene_symbol = data[1].strip()
                current_disease = data[2].strip().replace('"', "")
                new_disease = data[3].strip().replace('"', "")

                # Check if column "Updated" is defined in the file
                # If not, set value to the default "No"
                is_updated = data[4] if len(data) > 4 else "No"

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
                        print(
                            f"\n({g2p_id}) {gene_symbol}; {current_disease} (New disease: {new_disease})"
                        )

                    if current_disease in unique_diseases_from_input:
                        print(
                            f"WARNING: trying to update the same disease again: {current_disease}"
                        )
                    else:
                        unique_diseases_from_input.add(current_disease)

                    # Number of records linked to the gene symbol
                    n_records = 0

                    for record in db_data:
                        # In the old system most of the disease names don't include the '<gene>-related' in the name
                        # to compare disease names we have to compare the end of the name
                        if (
                            record["disease_name"].lower() == current_disease.lower()
                            or record["disease_name"].endswith(
                                "related " + current_disease
                            )
                            or record["disease_name"].endswith(
                                "associated " + current_disease
                            )
                        ) and record["disease_name"] != new_disease:
                            n_records += 1
                            to_update = 1

                            # Check if new disease name is already in the db
                            if new_disease in diseases:
                                # If the new disease name already exists then add it to replace existing name
                                # Before running the update, the API checks if records are not the same
                                print(
                                    f"Update LGD records: replace disease_id {record['disease_id']} by {diseases[new_disease]['disease_id']}"
                                )
                                lgd_disease_to_update.append(
                                    {
                                        "disease_id": record["disease_id"],
                                        "new_disease_id": diseases[new_disease][
                                            "disease_id"
                                        ],
                                    }
                                )

                                if dryrun:
                                    print(
                                        f"Update disease in lgd table -> replace disease_id: {record['disease_id']} by {diseases[new_disease]['disease_id']}"
                                    )

                            elif gene_symbol not in new_disease:
                                # This should not happen
                                # Action: print to file 'diseases_not_updated.txt'
                                wr_diseases.write(
                                    "Gene not found in new disease name\t" + line
                                )
                            else:
                                # Create list of diseases to update
                                diseases_to_update.append(
                                    {"id": record["disease_id"], "name": new_disease}
                                )
                                if dryrun:
                                    print(
                                        f"Update disease name -> disease_id: {record['disease_id']}; current name: {current_disease}; new name: {new_disease}"
                                    )

                        else:
                            list_disease.append(record["disease_name"])

                    # Flag if the disease is linked to more than one gene
                    if n_records > 1:
                        print(f"INFO: disease linked to {n_records} records")

                    # Print to file genes that won't have disease updates
                    if not to_update:
                        wr.write(
                            f"{g2p_id}\t{gene_symbol}\t{current_disease}\t{new_disease}\t{list_disease}\n"
                        )

    return diseases_to_update


def update_diseases(
    diseases_to_update: list, api_username: str, api_password: str, api_url: str
) -> None:
    """
    Method to update the disease names.
    This method calls the following gene2phenotype API endpoints:
    update/diseases/ to update the disease name (simple update)
    lgd_disease_updates/ to update the LGD record with the new disease id - it replaces the
    previous disease id with the new id

    Args:
        diseases_to_update (list): list of diseases to update
        api_username (str): G2P API username (super user)
        api_password (str): G2P API password
        api_url (str): URL to the G2P API
    """
    disease_url = "update/diseases/"
    lgd_disease_url = "lgd_disease_updates/"
    login_url = "login/"

    data = {"username": api_username, "password": api_password}

    response = requests.post(api_url + login_url, json=data)
    if response.status_code == 200:
        try:
            response_update = requests.post(
                api_url + disease_url, json=diseases_to_update, cookies=response.cookies
            )
            if response_update.status_code == 200:
                response_json = response_update.json()
                print("Diseases updated successfully:", response_json)
                # Some diseases were probably not updated: if they are already in the db
                # For those we can update the disease_id in lgd to point to the existing disease
                if "errors" in response_json:
                    for error in response_json["errors"]:
                        print(
                            f"Update LGD records: replace disease_id {error['id']} by {error['existing_id']}"
                        )
                        lgd_disease_to_update.append(
                            {
                                "disease_id": error["id"],
                                "new_disease_id": error["existing_id"],
                            }
                        )
            else:
                print(
                    "Failed to update diseases:",
                    response_update.status_code,
                    response_update.json(),
                )
        except Exception as e:
            print("Error while updating diseases:", e)

        if lgd_disease_to_update:
            try:
                response_update_lgd = requests.post(
                            api_url + lgd_disease_url,
                            json=lgd_disease_to_update,
                            cookies=response.cookies,
                        )
                if response_update_lgd.status_code == 200:
                    print(
                        "LGD records updated successfully:",
                        response_update_lgd.json(),
                    )
                else:
                    print(
                        "Failed to update LGD records:",
                        response_update_lgd.status_code,
                        response_update_lgd.json(),
                    )
            except Exception as e:
                print("Error while updating LGD diseases:", e)

    else:
        print("Error: cannot login into G2P")


def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--config", required=True, help="Config file")
    parser.add_argument(
        "--file",
        required=True,
        help="Tab delimited file with all diseases to be updated",
    )
    parser.add_argument(
        "--api_username", required=True, help="Username to connect to the G2P API"
    )
    parser.add_argument(
        "--api_password", required=True, help="Password to connect to the G2P API"
    )
    parser.add_argument(
        "--dryrun", required=False, default=0, help="Option to test the update"
    )
    args = parser.parse_args()

    file = args.file
    config_file = args.config
    dryrun = args.dryrun
    api_username = args.api_username
    api_password = args.api_password

    # Load the config file
    config = configparser.ConfigParser()
    config.read(config_file)

    db_host = config["database"]["host"]
    db_port = config["database"]["port"]
    db_name = config["database"]["name"]
    user = config["database"]["user"]
    password = config["database"]["password"]
    api_url = config["api"]["api_url"]

    print("Dump data from G2P...")
    gene_records, diseases = dump_data(db_host, int(db_port), db_name, user, password)
    print("Dump data from G2P... done\n")

    if os.path.isfile(file):
        expected_columns = [
            "g2p id",
            "gene symbol",
            "disease name",
            "disease name formatted",
            "allelic requirement", # NEW
        ]

        with open(file, "r", encoding="utf-8") as fh:
            header = fh.readline().strip().split("\t")
            if header[:4] != expected_columns:
                sys.exit(
                    f"Error: File format is incorrect. Found: {header[:4]}; Expected: {expected_columns}"
                )

        print("Parsing diseases to update...")
        diseases_to_update = read_file(file, gene_records, diseases, dryrun)
        print("Parsing diseases to update... done\n")

        if not dryrun:
            print("Updating disease names...")
            # This is a simple update of the disease name
            # The disease id stays the same - locus_genotype_disease still points to the same disease_id
            update_diseases(diseases_to_update, api_username, api_password, api_url)
            print("Updating disease names... done\n")

    else:
        print(f"Input file is invalid '{file}'")


if __name__ == "__main__":
    main()
