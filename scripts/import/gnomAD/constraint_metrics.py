import sys
import argparse
import MySQLdb
import datetime
import configparser

"""
    Script to import gene constraints metrics data from gnomAD into the gene_stats table of a G2P database.

    Params:
            --config : Config file name containing the database connection info (mandatory)
                    File format is the following:
                        [database]
                        host = <>
                        port = <>
                        user = <>
                        password = <>
                        name = <>

            --file : input file
                    Input file can be download from the URL: https://gnomad.broadinstitute.org/downloads#v4-constraint
                    File should be in TSV (.tsv) file format.

"""


# Constants
GNOMAD_SOURCE_NAME = "gnomAD"
ENSEMBL_SOURCE_NAME = "Ensembl"
PLI_ATTRIB_VALUE = "pli_gnomAD"
LOEUF_ATTRIB_VALUE = "loeuf_gnomAD"
META_KEY = "gnomAD_constraint_metrics"
META_DESCRIPTION = "gnomAD constraint metrics"
META_VERSION = "v4.1.0"


def get_attrib_id(
    attrib: str, db_host: str, db_port: int, db_name: str, user: str, password: str
) -> int:
    """
    Retrieve the ID corresponding to a specific attribute value from the 'attrib' table in the database.
    """
    get_attrib_id_query = """ SELECT id FROM attrib WHERE value = %s """

    database = MySQLdb.connect(
        host=db_host, port=db_port, user=user, passwd=password, db=db_name
    )
    cursor = database.cursor()

    cursor.execute(get_attrib_id_query, (attrib,))

    attrib_id = cursor.fetchone()

    if attrib_id is None:
        raise ValueError("Attrib ID not found in the database.")

    cursor.close()
    database.close()

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

    get_source_query = """ SELECT id FROM source WHERE name = %s """

    database = MySQLdb.connect(
        host=db_host, port=db_port, user=user, passwd=password, db=db_name
    )
    cursor = database.cursor()

    cursor.execute(get_source_query, (source_name,))
    source_id = cursor.fetchone()

    if source_id is None:
        raise ValueError(f"Source '{source_name}' not found in the database.")

    cursor.close()
    database.close()

    return source_id[0]


def validate_input_file(
    file: str,
) -> None:
    """
    Validate input file.

    Args:
        file (str): input file name
    Returns:
        None
    """
    try:
        with open(file, "r", encoding="utf-8") as fh:
            header = fh.readline().strip().split("\t")
            if header[1] != "gene_id":
                sys.exit(
                    f"Error: File format is incorrect. Found: {header[1]}; Expected: 'gene_id'"
                )
            elif header[3] != "canonical":
                sys.exit(
                    f"Error: File format is incorrect. Found: {header[3]}; Expected: 'canonical'"
                )
            elif header[18] != "lof.pLI":
                sys.exit(
                    f"Error: File format is incorrect. Found: {header[18]}; Expected: 'lof.pLI'"
                )
            elif header[22] != "lof.oe_ci.upper":
                sys.exit(
                    f"Error: File format is incorrect. Found: {header[22]}; Expected: 'lof.oe_ci.upper'"
                )
    except FileNotFoundError:
        sys.exit(f"File not found: {file}")
    except Exception as e:
        sys.exit(f"Error reading the file: {e}")


def is_valid_score(score: str) -> bool:
    """
    Validate if the score is valid.

    Args:
        score (str): score

    Returns:
        bool: if the score is valid or not
    """
    return score and score != "" and score != "NA" 


def get_ensembl_data_from_g2p_db(
    db_host: str, db_port: int, db_name: str, user: str, password: str
) -> dict[str, tuple[str, str]]:
    """
    Retrieves Ensembl ID to locus mappings from the G2P database.

    Args:
        db_host (str): hostname
        db_port (int): port
        db_name (str): G2P database name
        user (str): username
        password (str): password

    Returns:
        dict[str, tuple[str, str]]: dict of Ensembl ID and associated locus ID, locus name
    """
    ensembl_to_locus_mapping = {}

    get_ensembl_to_locus_mapping_query = """ SELECT li.identifier, li.locus_id, l.name FROM locus_identifier li 
                                                LEFT JOIN locus l ON li.locus_id = l.id 
                                                LEFT JOIN source s ON li.source_id = s.id 
                                                WHERE s.name = %s """

    database = MySQLdb.connect(
        host=db_host, port=db_port, user=user, passwd=password, db=db_name
    )
    cursor = database.cursor()

    cursor.execute(get_ensembl_to_locus_mapping_query, (ENSEMBL_SOURCE_NAME,))
    data = cursor.fetchall()
    for row in data:
        ensembl_id, locus_id, locus = row
        if ensembl_id not in ensembl_to_locus_mapping:
            ensembl_to_locus_mapping[ensembl_id] = (locus_id, locus)

    cursor.close()
    database.close()

    return ensembl_to_locus_mapping


def read_input_file(
    file: str,
    ensembl_to_locus_mapping: dict[str, tuple[str, str]],
    db_host: str,
    db_port: int,
    db_name: str,
    user: str,
    password: str,
) -> list[tuple[str, str, str, str, str]]:
    """
    Read and process input file.

    Args:
        file (str): input file name
        ensembl_to_locus_mapping(dict[str, tuple[str, str]]): dict of Ensembl ID and associated locus ID, locus name
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
    source_id = get_source_id(GNOMAD_SOURCE_NAME, db_host, db_port, db_name, user, password)

    # Get attrib id
    pli_attrib_id = get_attrib_id(
        PLI_ATTRIB_VALUE, db_host, db_port, db_name, user, password
    )
    loeuf_attrib_id = get_attrib_id(
        LOEUF_ATTRIB_VALUE, db_host, db_port, db_name, user, password
    )

    # input_ensembl_to_scores_mapping is used to store ensembl id (key) to scores (value) mappings in input data
    # The key-value mapping is of the following structure -> <ensembl_id>: (<pli_score>, <loeuf_score>)
    # IF a mapping is -> <ensembl_id>: () THEN it means that the <ensembl_id> does not have a canonical transcript in the input data
    # IF a mapping is -> <ensembl_id>: (<pli_score>, <loeuf_score>) THEN it means that the <ensembl_id> has a canonical transcript in the input data and the scores are stored
    input_ensembl_to_scores_mapping = {}
    
    # Read input file and store input data
    with open(file, "r", encoding="utf-8") as f:
        for line in f:
            if not line.startswith("gene"):
                line_data = line.split("\t")
                gene = line_data[0]
                gene_id = line_data[1]
                canonical = line_data[3]
                pli_score = line_data[18]
                loeuf_score = line_data[22]
                # Only rows with Ensembl gene IDS will be processed and others will be ignored
                if gene_id.startswith("ENSG"):
                    if gene_id not in input_ensembl_to_scores_mapping:
                        if canonical == "true":
                            input_ensembl_to_scores_mapping[gene_id] = (pli_score, loeuf_score)
                        elif canonical == "false":
                            input_ensembl_to_scores_mapping[gene_id] = ()
                    else:
                        if canonical == "true":
                            if len(input_ensembl_to_scores_mapping[gene_id]) == 0:
                                input_ensembl_to_scores_mapping[gene_id] = (pli_score, loeuf_score)
                            else:
                                print(
                                    f"Warning: Gene '{gene}' with Gene ID '{gene_id}' has multiple canonical transcripts in input file data. Processed first match but skipped row with Gene '{gene}', Gene ID '{gene_id}' and canonical value '{canonical}'."
                                )
    
    # Process input data
    success_count = 0
    for ensembl_id, score_tuple in input_ensembl_to_scores_mapping.items():
        if ensembl_id in ensembl_to_locus_mapping:
            if len(score_tuple) == 0:
                print(
                    f"Warning: Gene ID '{ensembl_id}' does not have canonical transcript in input file data. Skipped gene."
                )
            else:
                pli_score, loeuf_score = score_tuple
                locus_id, locus = ensembl_to_locus_mapping[ensembl_id]
                if is_valid_score(pli_score) and is_valid_score(loeuf_score):
                    final_data_to_insert.append(
                        (
                            locus,
                            locus_id,
                            pli_score,
                            source_id,
                            pli_attrib_id,
                        )
                    )
                    final_data_to_insert.append(
                        (
                            locus,
                            locus_id,
                            loeuf_score,
                            source_id,
                            loeuf_attrib_id,
                        )
                    )
                    success_count += 1
                elif is_valid_score(pli_score):
                    final_data_to_insert.append(
                        (
                            locus,
                            locus_id,
                            pli_score,
                            source_id,
                            pli_attrib_id,
                        )
                    )
                    print(
                        f"Warning: Gene ID '{ensembl_id}' has empty 'loeuf' score in input file data. Processed 'pli' score but skipped 'loeuf' score for this gene."
                    )
                    success_count += 1
                elif is_valid_score(loeuf_score):
                    final_data_to_insert.append(
                        (
                            locus,
                            locus_id,
                            loeuf_score,
                            source_id,
                            loeuf_attrib_id,
                        )
                    )
                    print(
                        f"Warning: Gene ID '{ensembl_id}' has empty 'pli' score in input file data. Processed 'loeuf' score but skipped 'pli' score for this gene."
                    )
                    success_count += 1
                else:
                    print(
                        f"Warning: Gene ID '{ensembl_id}' has empty 'pli' and 'loeuf' scores in input file data. Skipped gene."
                    )
        else:
            print(
                f"Warning: Gene ID '{ensembl_id}' not found in G2P DB. Skipped gene."
            )

    if len(final_data_to_insert) == 0:
        sys.exit("Error: No valid data found in input file.")

    # Print stats of input data processed
    print("----")
    print("Unique Ensembl Gene IDs found : ", len(input_ensembl_to_scores_mapping))
    print("Ensembl Gene IDs with valid data : ", success_count)
    print("Data rows to be inserted : ", len(final_data_to_insert))
    print("----")

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
    source_id = get_source_id(GNOMAD_SOURCE_NAME, db_host, db_port, db_name, user, password)

    insert_into_meta_query = """ INSERT into meta(`key`, date_update, description, version, source_id, is_public) VALUES (%s, %s, %s, %s, %s, %s) """

    database = MySQLdb.connect(
        host=db_host, port=db_port, user=user, passwd=password, db=db_name
    )
    cursor = database.cursor()

    cursor.execute(
        insert_into_meta_query,
        (META_KEY, datetime.datetime.now(), META_DESCRIPTION, META_VERSION, source_id, 0),
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
        help="Input file",
    )
    args = parser.parse_args()

    file = args.file
    config_file = args.config

    # Load the config file
    config = configparser.ConfigParser()
    config.read(config_file)

    db_host = config["database"]["host"]
    db_port = int(config["database"]["port"])
    db_name = config["database"]["name"]
    user = config["database"]["user"]
    pwd = config["database"]["password"]

    print("Validating input file...")
    validate_input_file(file)
    print("Validating input file... done\n")

    print("Getting ensembl ID to locus mapping from G2P DB...")
    ensembl_to_locus_mapping = get_ensembl_data_from_g2p_db(
        db_host, db_port, db_name, user, pwd
    )
    print("Getting ensembl ID to locus mapping from G2P DB... done\n")

    print("Reading input file...")
    final_scores = read_input_file(
        file,
        ensembl_to_locus_mapping,
        db_host,
        db_port,
        db_name,
        user,
        pwd,
    )
    print("Reading input file... done\n")

    print("Inserting data into gene_stats table...")
    insert_into_gene_stats(final_scores, db_host, db_port, db_name, user, pwd)
    print("Inserting data into gene_stats table... done\n")

    print("Inserting data into meta table...")
    insert_details_into_meta(db_host, db_port, db_name, user, pwd)
    print("Inserting data into meta table... done\n")


if __name__ == "__main__":
    main()
