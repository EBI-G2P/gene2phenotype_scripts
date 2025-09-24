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
                    g2p id\tgene symbol\tdisease name\tdisease name formatted\tallelic requirement\tUpdated

            --api_username: Username to connect to the G2P API (mandatory)
            --api_password: Password to connect to the G2P API (mandatory)
            --dryrun: Test script without running the updates (not run by default)
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
        dict[str, dict]: dict of diseases by name and associated records
    """
    diseases = {}  # key = disease name; value = list of record ids

    sql_disease = """ SELECT d.name, d.id, lgd.id, g.stable_id, g.is_deleted, a.value
                      FROM disease d
                      LEFT JOIN locus_genotype_disease lgd ON lgd.disease_id = d.id
                      LEFT JOIN g2p_stableid g ON g.id = lgd.stable_id
                      LEFT JOIN attrib a ON a.id = lgd.genotype_id
                      ORDER BY d.name
                  """

    db = MySQLdb.connect(
        host=db_host, port=db_port, user=user, passwd=password, db=db_name
    )
    cursor = db.cursor()
    cursor.execute(sql_disease)
    data_disease = cursor.fetchall()
    for row in data_disease:
        if row[0] not in diseases:
            diseases[row[0]] = {
                "disease_id": row[1],
                "records": [
                    {
                        "lgd_id": row[2],
                        "stable_id": row[3],
                        "record_is_deleted": row[4],
                        "genotype": row[5],
                    }
                ],
            }
        else:
            diseases[row[0]]["records"].append(
                {
                    "lgd_id": row[2],
                    "stable_id": row[3],
                    "record_is_deleted": row[4],
                    "genotype": row[5],
                }
            )

    db.close()

    return diseases


def login(
    api_username: str, api_password: str, api_url: str
) -> requests.cookies.RequestsCookieJar:
    """Login into G2P API"""
    login_url = f"{api_url.rstrip('/')}/login/"

    response = requests.post(
        login_url, json={"username": api_username, "password": api_password}
    )

    if response.status_code != 200:
        sys.exit("Login failed. Check your credentials and API URL.")

    return response.cookies


def logout(api_url: str, cookies: requests.cookies.RequestsCookieJar) -> None:
    """Logout of the API"""
    logout_url = f"{api_url.rstrip('/')}/logout/"

    response = requests.post(logout_url, cookies=cookies)

    if response.status_code != 204:
        sys.exit("Logout failed. Check your credentials and API URL.")


def read_file(
    file: str,
    diseases: dict[str, dict],
    dryrun: bool,
    api_url: str,
    cookies: requests.cookies.RequestsCookieJar,
) -> list:
    """
    Reads the diseases from the input file, for each disease checks if it can be udpdated.

    Args:
        file (str): tab delimited input file
        diseases (dict[str, dict]): dictionary of diseases by name and associated records
        dryrun (bool): run script in update mode (default: 0)
        api_url (str): G2P API URL
        cookies (requests.cookies.RequestsCookieJar): cookies contained tokens (login credentials)

    Returns:
        list: list of diseases to be updated
    """
    file_not_updated = "diseases_not_updated.txt"
    diseases_to_update = []

    with (
        open(file, "r") as fh,
        open(file_not_updated, "w") as wr_diseases,
    ):
        wr_diseases.write(
            "comment\tg2p id\tgene symbol\tdisease name\tdisease name formatted\tgenotype\n"
        )

        # Header:
        # g2p id, gene symbol, disease name, disease name formatted, allelic requirement, Updated
        # Other columns are ignored
        for line in fh:
            if not line.startswith("g2p id"):
                data = line.strip().split("\t")
                g2p_id = data[0].strip()
                gene_symbol = data[1].strip()
                current_disease = data[2].strip().replace('"', "")
                new_disease = data[3].strip().replace('"', "")
                genotype = data[4].strip()

                # Check if column "Updated" is defined in the file
                # If not, set value to the default "No"
                is_updated = data[5] if len(data) > 5 else "No"

                if is_updated == "Yes":
                    continue

                # Check if the record can be found in G2P
                try:
                    db_data = diseases[current_disease]
                except KeyError:
                    wr_diseases.write(
                        f"Current disease not found in G2P\t{g2p_id}\t{gene_symbol}\t{current_disease}\t{new_disease}\t{genotype}\n"
                    )
                else:
                    to_update = 0

                    if dryrun:
                        print(
                            f"\n{g2p_id}; {gene_symbol}; {genotype}; {current_disease} (New disease: {new_disease})"
                        )

                    # Check if new disease follows dyadic name with correct gene symbol
                    if not new_disease.startswith(gene_symbol + "-related"):
                        wr_diseases.write(
                            f"Invalid new disease name\t{g2p_id}\t{gene_symbol}\t{current_disease}\t{new_disease}\t{genotype}\n"
                        )
                        continue

                    # Check if we are trying to update the same disease again
                    # Note: use the disease and genotype - the same disease can be updated for different genotypes
                    if current_disease + "-" + genotype in unique_diseases_from_input:
                        wr_diseases.write(
                            f"This disease + genotype has already been updated\t{g2p_id}\t{gene_symbol}\t{current_disease}\t{new_disease}\t{genotype}\n"
                        )
                        continue
                    else:
                        unique_diseases_from_input.add(current_disease + "-" + genotype)

                    # Number of records linked to the current disease
                    n_records = len(db_data["records"])

                    if n_records > 1:
                        # CASE 1: disease associated with several records
                        # Get from the db data the record that has the same genotype
                        found = 0
                        for record_from_db in db_data["records"]:
                            if (
                                record_from_db["genotype"] == genotype
                                and not record_from_db["record_is_deleted"]
                            ):
                                found = 1
                                # Check if new disease name exists in G2P
                                try:
                                    db_data_new_disease = diseases[new_disease]
                                except KeyError:
                                    # If disease name is not in G2P then call endpoint to create disease
                                    # Then update the LGD record to point to new disease name
                                    if dryrun:
                                        print(f"   Going to add disease: {new_disease}")
                                        print(
                                            f"   Update disease in lgd table. Going to update LGD record {g2p_id}: replace disease_id {db_data['disease_id']} by new disease {new_disease}"
                                        )
                                    else:
                                        # Call endpoint to create new disease
                                        new_disease_id = add_disease(
                                            new_disease, api_url, cookies
                                        )
                                        if not new_disease_id:
                                            wr_diseases.write(
                                                f"Could not insert new disease name\t{g2p_id}\t{gene_symbol}\t{current_disease}\t{new_disease}\t{genotype}\n"
                                            )
                                            continue
                                        lgd_disease_to_update.append(
                                            {
                                                "disease_id": db_data["disease_id"],
                                                "new_disease_id": new_disease_id,
                                                "stable_id": g2p_id,
                                            }
                                        )
                                else:
                                    if dryrun:
                                        print(
                                            f"   Update disease in lgd table: disease '{new_disease}' already exists in G2P. Going to update LGD record {g2p_id}: replace disease_id {db_data['disease_id']} by {db_data_new_disease['disease_id']}"
                                        )
                                    lgd_disease_to_update.append(
                                        {
                                            "disease_id": db_data["disease_id"],
                                            "new_disease_id": db_data_new_disease[
                                                "disease_id"
                                            ],
                                            "stable_id": g2p_id,
                                        }
                                    )
                        if not found:
                            wr_diseases.write(
                                f"Record cannot be found in G2P (maybe it's deleted)\t{g2p_id}\t{gene_symbol}\t{current_disease}\t{new_disease}\t{genotype}\n"
                            )
                    else:
                        # CASE 2: Disease is associated with only 1 record - we can update it only if record not deleted
                        if db_data["records"][0]["record_is_deleted"]:
                            print(f"WARNING: {g2p_id} is deleted")
                            wr_diseases.write(
                                f"{g2p_id} is deleted\t{g2p_id}\t{gene_symbol}\t{current_disease}\t{new_disease}\t{genotype}\n"
                            )
                            continue

                        # Check if the record id is the same as the g2p id from the input file
                        if g2p_id != db_data["records"][0]["stable_id"]:
                            if dryrun:
                                print(
                                    f"WARNING: {g2p_id} looks suspicious. Please check record {db_data['records'][0]['stable_id']}"
                                )
                            wr_diseases.write(
                                f"{g2p_id} looks suspicious. Please check record {db_data['records'][0]['stable_id']}\t{g2p_id}\t{gene_symbol}\t{current_disease}\t{new_disease}\t{genotype}\n"
                            )
                            continue

                        try:
                            db_data_new_disease = diseases[new_disease]
                        except KeyError:
                            to_update = 1
                        else:
                            # New disease name is already in the db
                            # Add disease to list of lgd to update
                            if dryrun:
                                print(
                                    f"Update disease in lgd table: disease '{new_disease}' already exists in G2P. Going to update LGD record {g2p_id}: replace disease_id {db_data['disease_id']} by {db_data_new_disease['disease_id']}"
                                )
                            lgd_disease_to_update.append(
                                {
                                    "disease_id": db_data["disease_id"],
                                    "new_disease_id": db_data_new_disease["disease_id"],
                                    "stable_id": g2p_id,
                                }
                            )

                        # Add disease to list of diseases to update
                        if to_update:
                            diseases_to_update.append(
                                {"id": db_data["disease_id"], "name": new_disease}
                            )
                            if dryrun:
                                print(
                                    f"Update disease name -> disease_id: {db_data['disease_id']}; current name: {current_disease}; new name: {new_disease}"
                                )

    return diseases_to_update


def add_disease(
    disease_name: str, api_url: str, cookies: requests.cookies.RequestsCookieJar
) -> int:
    """
    Method to create a new disease.
    This method calls the following gene2phenotype API endpoint:
        - add/disease/
    It returns the disease id just that was just created.
    """
    add_disease_url = "add/disease/"
    disease_id = None

    # Prepare disease input data
    disease_data = {"name": disease_name}

    try:
        response_add = requests.post(
            api_url + add_disease_url, json=disease_data, cookies=cookies
        )
        if response_add.status_code == 201:
            response_json = response_add.json()
            disease_id = response_json["id"]
            print(f"Disease added -> id: {disease_id}; {response_json['name']}")
        else:
            print(
                "Failed to add disease:",
                disease_name,
                response_add.json(),
            )
    except Exception as e:
        print("Error while adding disease:", str(e))

    return disease_id


def update_diseases(
    diseases_to_update: list, api_url: str, cookies: requests.cookies.RequestsCookieJar
) -> None:
    """
    Method to update the disease names.
    This method calls the following gene2phenotype API endpoints:
    update/diseases/ to update the disease name (simple update)
    lgd_disease_updates/ to update the LGD record with the new disease id - it replaces the
    previous disease id with the new id

    Args:
        diseases_to_update (list): list of diseases to update
    """
    disease_url = "update/diseases/"
    lgd_disease_url = "lgd_disease_updates/"

    try:
        response_update = requests.post(
            api_url + disease_url, json=diseases_to_update, cookies=cookies
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
        print("Error while updating diseases:", str(e))

    if lgd_disease_to_update:
        try:
            response_update_lgd = requests.post(
                api_url + lgd_disease_url,
                json=lgd_disease_to_update,
                cookies=cookies,
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
        "--dryrun",
        action="store_true",
        help="Option to test the which diseases are going to be updated",
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
    cookies = None

    print("Dump data from G2P...")
    diseases = dump_data(db_host, int(db_port), db_name, user, password)
    print("Dump data from G2P... done\n")

    if os.path.isfile(file):
        expected_columns = [
            "g2p id",
            "gene symbol",
            "disease name",
            "disease name formatted",
            "allelic requirement",
        ]

        with open(file, "r", encoding="utf-8") as fh:
            header = fh.readline().strip().split("\t")
            if header[:5] != expected_columns:
                sys.exit(
                    f"Error: File format is incorrect. Found: {header[:5]}; Expected: {expected_columns}"
                )

        if not dryrun:
            print("Logging in...")
            cookies = login(api_username, api_password, api_url)

        print("Parsing diseases to update...")
        diseases_to_update = read_file(file, diseases, dryrun, api_url, cookies)
        print("Parsing diseases to update... done\n")

        if not dryrun:
            print("Updating disease names...")
            update_diseases(diseases_to_update, api_url, cookies)
            print("Updating disease names... done\n")

            print("Logging out...")
            logout(api_url, cookies)
        else:
            print("Dry run: disease updates skipped")

    else:
        print(f"Input file is invalid '{file}'")


if __name__ == "__main__":
    main()
