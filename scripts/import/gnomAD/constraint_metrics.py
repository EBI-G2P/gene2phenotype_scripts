import os.path
import sys
import argparse
import MySQLdb
import datetime
import configparser

"""
    Script to import gene constraints metrics data from gnomAD into the gene_stats table of a G2P database.

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

            --file : input file path

"""


# Constants
SOURCE_NAME = "gnomAD constraint metrics"
PLI_ATTRIB_VALUE = "pli_gnomAD"
LOEUF_ATTRIB_VALUE = "loeuf_gnomAD"


def get_locus_id_from_g2p_db(
    db_host: str, db_port: int, db_name: str, user: str, password: str
) -> dict[str, str]:
    """
    Retrieves locus IDs from the G2P database.

    Args:
        db_host (str): hostname
        db_port (int): port
        db_name (str): G2P database name
        user (str): username
        password (str): password

    Returns:
        dict[str, str]: dict of genes and associated ids
    """
    final_list = {}

    sql_get_locus = """ SELECT id, name from locus """
    sql_get_locus_attrib = """ SELECT locus_id, value from locus_attrib """

    database = MySQLdb.connect(
        host=db_host, port=db_port, user=user, passwd=password, db=db_name
    )
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


def read_input_file(
    file: str,
    list_gene_ids: dict[str, str],
    db_host: str,
    db_port: int,
    db_name: str,
    user: str,
    password: str,
) -> list[tuple[str, str, str, str, str]]:
    """
    Read input file.

    Args:
        file (str): input file name
        list_gene_ids(dict[str, str]): dict of genes and associated ids
        db_host (str): hostname
        db_port (int): port
        db_name (str): G2P database name
        user (str): username
        password (str): password

    Returns:
        list[tuple[str, str, str, str, str]]: list of tuples containing db values to insert
    """
    final_data_to_insert = []

    # Get source id
    source_id = get_source_id(SOURCE_NAME, db_host, db_port, db_name, user, password)

    # Get attrib id
    pli_attrib_id = get_attrib_id(
        PLI_ATTRIB_VALUE, db_host, db_port, db_name, user, password
    )
    loeuf_attrib_id = get_attrib_id(
        LOEUF_ATTRIB_VALUE, db_host, db_port, db_name, user, password
    )

    f = open(file, "r")
    previous_gene = ""
    is_previous_gene_mane_select_transcript_exists = False
    for line in f:
        if not line.startswith("gene"):
            line_data = line.split("\t")
            gene = line_data[0]
            if previous_gene != "" and gene != previous_gene:
                if not is_previous_gene_mane_select_transcript_exists:
                    print(f"Gene '{gene}' does not have MANE transcript.")
                previous_gene = gene
                is_previous_gene_mane_select_transcript_exists = False
            if gene in list_gene_ids:
                mane_select = line_data[4]
                if mane_select == "TRUE":
                    pli_score = line_data[18]
                    pli_score_line = (
                        gene,
                        list_gene_ids[gene],
                        pli_score,
                        source_id,
                        pli_attrib_id,
                    )
                    final_data_to_insert.append(pli_score_line)
                    loeuf_score = line_data[22]
                    loeuf_score_line = (
                        gene,
                        list_gene_ids[gene],
                        loeuf_score,
                        source_id,
                        loeuf_attrib_id,
                    )
                    final_data_to_insert.append(loeuf_score_line)
                    is_previous_gene_mane_select_transcript_exists = True
            else:
                print(f"Gene '{gene}' not found in G2P DB.")
    if previous_gene != "":
        if not is_previous_gene_mane_select_transcript_exists:
            print(f"Gene '{gene}' does not have MANE transcript.")
    f.close()

    return final_data_to_insert


def insert_into_gene_stats(
    final_scores: list[tuple[str, str, str, str, str]],
    db_host: str,
    db_port: int,
    db_name: str,
    user: str,
    password: str,
) -> None:
    """
    Inserts gene-related data into the `gene_stats` table.

    Args:
        final_scores (list[tuple[str, str, str, str, str]]): list of tuples containing db values to insert
        db_host (str): hostname
        db_port (int): port
        db_name (str): G2P database name
        user (str): username
        password (str): password

    Returns:
        None
    """

    sql_insert_gene_stats = """ INSERT into gene_stats (gene_symbol, gene_id, score, source_id, description_attrib_id) VALUES (%s, %s, %s, %s, %s)"""

    database = MySQLdb.connect(
        host=db_host, port=db_port, user=user, passwd=password, db=db_name
    )
    cursor = database.cursor()
    cursor.executemany(sql_insert_gene_stats, final_scores)
    database.commit()
    cursor.close()
    database.close()


def get_attrib_id(
    attrib: str, db_host: str, db_port: int, db_name: str, user: str, password: str
) -> int:
    """
    Retrieve the ID corresponding to a specific attribute value from the 'attrib' table in the database.
    """
    get_attrib_id_query = """SELECT id from attrib where value = %s"""

    database = MySQLdb.connect(
        host=db_host, port=db_port, user=user, passwd=password, db=db_name
    )
    cursor = database.cursor()

    cursor.execute(get_attrib_id_query, (attrib,))

    attrib_id = cursor.fetchone()

    if attrib_id is None:
        raise ValueError("Attrib ID not found in the database.")

    return attrib_id[0]


def get_source_id(
    source_name: str, db_host: str, db_port: int, db_name: str, user: str, password: str
) -> int:
    """
    Retrieves the source ID corresponding to a specific source value from the 'source' table in the database.

    Args:
        source_name (str): source name
        db_host (str): hostname
        db_port (int): port
        db_name (str): G2P database name
        user (str): username
        password (str): password

    Returns:
        int: source ID
    """

    get_source_query = """ SELECT id from source where name = %s """

    database = MySQLdb.connect(
        host=db_host, port=db_port, user=user, passwd=password, db=db_name
    )
    cursor = database.cursor()

    cursor.execute(get_source_query, (source_name,))
    source_id = cursor.fetchone()

    if source_id is None:
        raise ValueError("Source ID not found in the database.")

    source_id = source_id[0]

    cursor.close()
    database.close()

    return source_id


def insert_details_into_meta(
    db_host: str, db_port: int, db_name: str, user: str, password: str
) -> None:
    """
    Inserts data import details into the 'meta' table.

    Args:
        db_host (str): hostname
        db_port (int): port
        db_name (str): G2P database name
        user (str): username
        password (str): password

    Returns:
        None
    """
    source_id = get_source_id(SOURCE_NAME, db_host, db_port, db_name, password, user)

    insert_into_meta_query = """ INSERT into meta(key, date_update, description, version, source_id, is_public) VALUES (%s, %s, %s, %s, %s, %s) """

    database = MySQLdb.connect(
        host=db_host, port=db_port, user=user, passwd=password, db=db_name
    )
    cursor = database.cursor()

    meta_key = "gnomAD_constraint_metrics"
    description = "gnomAD constraint metrics"
    version = "v4"
    cursor.execute(
        insert_into_meta_query,
        (meta_key, datetime.datetime.now(), description, version, source_id, 0),
    )

    database.commit()
    cursor.close()
    database.close()


def main():
    parser = argparse.ArgumentParser(
        description="This script is used to import the constraint metrics from a file and insert that data into the 'gene_stats' table in the G2P DB"
    )
    parser.add_argument("--config", required=True, help="Config file")
    parser.add_argument(
        "--file",
        required=True,
        help="Input file path",
    )
    args = parser.parse_args()

    file = args.file
    config_file = args.config

    # Load the config file
    config = configparser.ConfigParser()
    config.read(config_file)

    db_host = config["database"]["host"]
    db_port = config["database"]["port"]
    db_name = config["database"]["name"]
    user = config["database"]["user"]
    pwd = config["database"]["password"]

    if os.path.isfile(file):
        # Validate input file
        with open(file, "r", encoding="utf-8") as fh:
            header = fh.readline().strip().split("\t")
            if header[0] != "gene":
                sys.exit(
                    f"Error: File format is incorrect. Found: {header[0]}; Expected: 'gene'"
                )
            elif header[4] != "mane_select":
                sys.exit(
                    f"Error: File format is incorrect. Found: {header[4]}; Expected: 'mane_select'"
                )
            elif header[18] != "lof.pLI":
                sys.exit(
                    f"Error: File format is incorrect. Found: {header[18]}; Expected: 'lof.pLI'"
                )
            elif header[22] != "lof.oe_ci.upper":
                sys.exit(
                    f"Error: File format is incorrect. Found: {header[22]}; Expected: 'lof.oe_ci.upper'"
                )

        print("Getting locus id from G2P DB...")
        list_gene_ids = get_locus_id_from_g2p_db(
            db_host, int(db_port), db_name, user, pwd
        )
        print("Getting locus id from G2P DB... done\n")

        print("Reading input file...")
        final_scores = read_input_file(
            file, list_gene_ids, db_host, db_port, db_name, user, pwd
        )
        print("Reading input file... done\n")

        print("Inserting data into gene_stats table...")
        insert_into_gene_stats(final_scores, db_host, db_port, db_name, user, pwd)
        print("Inserting data into gene_stats table... done\n")

        print("Inserting data into meta table...")
        insert_details_into_meta(db_host, db_port, db_name, user, pwd)
        print("Inserting data into meta table... done\n")
    else:
        print(f"Input file '{file}' is invalid")


if __name__ == "__main__":
    main()
