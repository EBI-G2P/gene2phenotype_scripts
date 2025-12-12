#!/usr/bin/env python3

import argparse
import configparser
import datetime
import os
import sys

import MySQLdb


"""
    Script to import gene probability data into the gene_stats table of a G2P database.

    Options:
            --config : Config file name containing the G2P database connection info (mandatory)
                            File format is the following:
                                [g2p_database]
                                host = <>
                                port = <>
                                user = <>
                                password = <>
                                name = <>

            -d, --dir : str (mandatory)
                Path to the folder with all the input files.

    Example usage:
            python script.py --config config.ini --dir /path/to/files
"""


# Mapping the attrib to the input file
attrib_mapping = {
    "gain_of_function.tsv": "gain_of_function_mp",
    "loss_of_function.tsv": "loss_of_function_mp",
    "dominant_negative.tsv": "dominant_negative_mp",
}

KEY = "Badonyi_score"


def get_locus_id_from_g2p_db(host, port, db, password, user):
    """
    Retrieves locus IDs from the G2P database.
    """
    final_list = {}

    sql_get_locus = """ SELECT id, name from locus """
    sql_get_locus_attrib = """SELECT locus_id, value from locus_attrib """

    database = MySQLdb.connect(host=host, port=port, user=user, passwd=password, db=db)
    cursor = database.cursor()

    cursor.execute(sql_get_locus)
    data = cursor.fetchall()
    for row in data:
        if row[1] not in final_list:
            final_list[row[1]] = row[0]

    cursor.execute(sql_get_locus_attrib)
    data_2 = cursor.fetchall()
    for row in data_2:
        if row[1] not in final_list:
            final_list[row[1]] = row[0]

    cursor.close()
    database.close()

    return final_list


def read_input_files(dir, list_gene_ids, host, port, db, password, user):
    final_data_to_insert = []

    source_id = get_source_details(host, port, db, password, user)

    files = os.listdir(dir)
    for file in files:
        if os.path.isfile(os.path.join(dir, file)):
            # Get the attrib value based on the file name
            attrib_value = attrib_mapping[file]
            attrib_id = get_attrib_ids(host, port, db, password, user, attrib_value)

            f = open(os.path.join(dir, file), "r")
            for line in f:
                if not line.startswith("gene"):
                    line_data = line.split("\t")
                    # check if gene name is in G2P db
                    if line_data[0] in list_gene_ids:
                        final_line = (
                            line_data[0],
                            list_gene_ids[line_data[0]],
                            line_data[2],
                            source_id,
                            attrib_id,
                        )
                        final_data_to_insert.append(final_line)

            f.close()

    return final_data_to_insert


def insert_into_gene_stats(host, port, db, password, user, final_scores):
    """
    Reads input files and inserts gene-related data into the `gene_stats` table.
    """

    sql_insert_gene_stats = """ INSERT into gene_stats (gene_symbol, gene_id, score, source_id, description_attrib_id)
                                VALUES (%s, %s, %s, %s, %s)"""

    database = MySQLdb.connect(host=host, port=port, user=user, passwd=password, db=db)
    cursor = database.cursor()
    cursor.executemany(sql_insert_gene_stats, final_scores)
    database.commit()
    cursor.close()
    database.close()


def get_attrib_ids(host, port, db, password, user, attrib):
    """
    Retrieve the ID corresponding to a specific attribute value from a MySQL database.
    """
    get_attrib_id_query = """SELECT id from attrib where value = %s"""

    database = MySQLdb.connect(host=host, port=port, user=user, passwd=password, db=db)
    cursor = database.cursor()

    cursor.execute(get_attrib_id_query, (attrib,))

    attrib_id = cursor.fetchone()

    if attrib_id is None:
        raise ValueError("Source ID not found in the database.")

    return attrib_id[0]


def get_source_details(host, port, db, password, user):
    """
    Retrieves the source ID from the 'source' table in the database.
    """

    get_source_query = (
        """ SELECT id from source where name = 'Marsh Mechanism probabilities'"""
    )

    database = MySQLdb.connect(host=host, port=port, user=user, passwd=password, db=db)
    cursor = database.cursor()

    cursor.execute(get_source_query)
    source_id = cursor.fetchone()

    if source_id is None:
        raise ValueError("Source ID not found in the database.")

    source_id = source_id[0]

    cursor.close()
    database.close()

    return source_id


def insert_details_into_meta(host, port, db, password, user):
    """
    Inserts import details into the 'meta' table.
    """

    source_id = get_source_details(host, port, db, password, user)
    description = "Baydoni & Marsh probabilities"

    insert_into_meta_query = """ INSERT into meta(`key`, date_update, description, version, source_id, is_public)
                                 VALUES (%s, %s, %s, %s, %s, %s) """

    database = MySQLdb.connect(host=host, port=port, user=user, passwd=password, db=db)
    cursor = database.cursor()

    for file, attrib in attrib_mapping.items():
        # current_datetime = datetime.now()
        meta_key = KEY + "_" + attrib
        # formatted_datetime = current_datetime.strftime('%Y-%m-%d %H:%M:%S.%f')
        version = 1
        cursor.execute(
            insert_into_meta_query,
            (meta_key, datetime.datetime.now(), description, version, source_id, 0),
        )

    database.commit()
    cursor.close()
    database.close()


def main():
    parser = argparse.ArgumentParser(
        description="This script is used to import the probabilities from a file and imports it to the gene_stats table in the G2P DB"
    )
    parser.add_argument(
        "--config",
        required=True,
        help="Config file with G2P database connection details",
    )
    parser.add_argument("--dir", required=True, help="Directory with the input files")

    args = parser.parse_args()

    config_file = args.config
    dir = args.dir

    # Load the config file
    config = configparser.ConfigParser()
    config.read(config_file)

    try:
        g2p_config = config["g2p_database"]
    except KeyError:
        sys.exit("ERROR: 'g2p_database' missing from config file")
    else:
        host = g2p_config["host"]
        port = int(g2p_config["port"])
        db = g2p_config["name"]
        user = g2p_config["user"]
        pwd = g2p_config["password"]

    print("Getting locus id from the G2P DB...")
    list_gene_ids = get_locus_id_from_g2p_db(host, port, db, pwd, user)
    print("Getting locus id from the G2P DB... done\n")

    print("Reading input file...")
    final_scores = read_input_files(dir, list_gene_ids, host, port, db, pwd, user)
    print("Reading input file... done\n")

    print("Inserting into gene stats...")
    insert_into_gene_stats(host, port, db, pwd, user, final_scores)
    print("Inserting into gene stats... done\n")

    insert_details_into_meta(host, port, db, pwd, user)
    print("Meta updated")
    sys.exit


if __name__ == "__main__":
    main()
