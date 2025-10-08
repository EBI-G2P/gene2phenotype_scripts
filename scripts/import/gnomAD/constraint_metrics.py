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
            if header[0] != "gene":
                sys.exit(
                    f"Error: File format is incorrect. Found: {header[0]}; Expected: 'gene'"
                )
            elif header[1] != "gene_id":
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


def validate_previous_gene_data(
    locus_to_locus_id_mapping: dict[str, str],
    ensembl_to_gene_mapping: dict[str, str],
    previous_gene: str,
    previous_ensembl_gene_id: str,
    is_success: bool,
) -> None:
    """
    Validate previous gene data. This function is used when iterating over input file rows.

    Args:
        locus_to_locus_id_mapping(dict[str, str]): dict of genes and associated ids
        ensembl_to_gene_mapping(dict[str, str]): dict of ensembl id and associated gene name
        previous_gene (str): previous gene
        previous_ensembl_gene_id (str): previous Ensembl gene id
        is_success (bool): Flag to track if previous gene data was a success

    Returns:
        None
    """
    if previous_gene != "" and not is_success:
        if previous_gene not in locus_to_locus_id_mapping:
            if previous_ensembl_gene_id == "":
                print(
                    f"Warning: Gene '{previous_gene}' not found in G2P DB. Skipped gene."
                )
            elif previous_ensembl_gene_id not in ensembl_to_gene_mapping:
                print(
                    f"Warning: Gene '{previous_gene}' with Ensembl Gene Id '{previous_ensembl_gene_id}' not found in G2P DB. Skipped gene."
                )
            elif previous_ensembl_gene_id in ensembl_to_gene_mapping:
                print(
                    f"Warning: Gene '{previous_gene}' not found in G2P DB but Ensembl Gene Id '{previous_ensembl_gene_id}' found in G2P DB linked to different Gene '{ensembl_to_gene_mapping[previous_ensembl_gene_id]}'. Skipped gene."
                )
        elif previous_gene in locus_to_locus_id_mapping:
            if previous_ensembl_gene_id == "":
                print(
                    f"Warning: Gene '{previous_gene}' does not have associated Ensembl data in input file data. Skipped gene."
                )
            elif not is_success:
                print(
                    f"Warning: Gene '{previous_gene}' does not have canonical transcript in input file data. Skipped gene."
                )


def validate_scores(
    locus_to_locus_id_mapping: dict[str, str],
    gene: str,
    gene_id: str,
    pli_score: str,
    loeuf_score: str,
    pli_attrib_id: str,
    loeuf_attrib_id: str,
    source_id: str,
    is_empty_gene: bool = False,
) -> list[tuple[str, str, str, str, str]]:
    """
    Validate pli and loeuf scores. This function is used when iterating over input file rows.

    Args:
        locus_to_locus_id_mapping(dict[str, str]): dict of gene and associated ID
        gene (str): gene
        gene_id (str): gene ID
        pli_score (str): PLI score
        loeuf_score (str):  LOEUF score
        pli_attrib_id (str): PLI Attrib ID in G2P DB
        loeuf_attrib_id (str): LOEUF Attrib ID in G2P DB
        source_id (str): Source ID in G2P DB
        is_empty_gene (bool): Boolean flag to check if gene is 'NA' gene

    Returns:
        list[tuple[str, str, str, str, str]]: list of tuples containing db values to insert
    """
    score_data_to_insert = []
    gene_to_print = ""
    if is_empty_gene:
        gene_to_print = "NA"
    else:
        gene_to_print = gene
    if (not pli_score or pli_score == "NA" or pli_score == "") and (
        not loeuf_score or loeuf_score == "NA" or loeuf_score == ""
    ):
        print(
            f"Warning: Gene '{gene_to_print}' with Gene ID '{gene_id}' has empty 'pli' and 'loeuf' scores in input file data. Skipped gene."
        )
    elif not loeuf_score or loeuf_score == "NA" or loeuf_score == "":
        score_data_to_insert.append(
            (
                gene,
                locus_to_locus_id_mapping[gene],
                pli_score,
                source_id,
                pli_attrib_id,
            )
        )
        print(
            f"Warning: Gene '{gene_to_print}' has empty 'loeuf' score in input file data. Processed 'pli' score but skipped 'loeuf' score."
        )
    elif not pli_score or pli_score == "NA" or pli_score == "":
        score_data_to_insert.append(
            (
                gene,
                locus_to_locus_id_mapping[gene],
                loeuf_score,
                source_id,
                loeuf_attrib_id,
            )
        )
        print(
            f"Warning: Gene '{gene_to_print}' has empty 'pli' score in input file data. Processed 'loeuf' score but skipped 'pli' score."
        )
    else:
        score_data_to_insert.append(
            (
                gene,
                locus_to_locus_id_mapping[gene],
                pli_score,
                source_id,
                pli_attrib_id,
            )
        )
        score_data_to_insert.append(
            (
                gene,
                locus_to_locus_id_mapping[gene],
                loeuf_score,
                source_id,
                loeuf_attrib_id,
            )
        )
    return score_data_to_insert


def process_empty_genes(
    empty_genes_data_list: list[tuple[str, str, str, str]],
    ensembl_to_gene_mapping: dict[str, str],
    locus_to_locus_id_mapping: dict[str, str],
    pli_attrib_id: str,
    loeuf_attrib_id: str,
    source_id: str,
):
    """
    Process 'NA' genes data

    Args:
        empty_genes_data_list(list[tuple[str, str, str, str]]): data with empty genes data
        ensembl_to_gene_mapping(dict[str, str]): dict of ensembl id and associated gene name
        locus_to_locus_id_mapping(dict[str, str]): dict of gene and associated ID
        pli_attrib_id (str): PLI Attrib ID in G2P DB
        loeuf_attrib_id (str): LOEUF Attrib ID in G2P DB
        source_id (str): Source ID in G2P DB

    Returns:
        list[tuple[str, str, str, str, str]]: list of tuples containing db values to insert
    """
    data_to_insert = []
    ensembl_gene_id_to_canonical_mapping = {}
    non_ensembl_gene_id_list = []
    not_found_ensembl_gene_id_list = []
    # Process 'NA' Genes data
    for item in empty_genes_data_list:
        gene_id, canonical, pli_score, loeuf_score = item
        if not gene_id.startswith("ENSG"):
            if gene_id not in non_ensembl_gene_id_list:
                print(
                    f"Warning: Gene 'NA' with Gene ID '{gene_id}' is not Ensembl Gene Id. Skipped gene."
                )
                non_ensembl_gene_id_list.append(gene_id)
        else:
            if gene_id in ensembl_to_gene_mapping:
                if gene_id not in ensembl_gene_id_to_canonical_mapping:
                    if canonical == "true":
                        gene = ensembl_to_gene_mapping[gene_id]
                        score_data_to_insert = validate_scores(
                            locus_to_locus_id_mapping,
                            gene,
                            gene_id,
                            pli_score,
                            loeuf_score,
                            pli_attrib_id,
                            loeuf_attrib_id,
                            source_id,
                            True,
                        )
                        data_to_insert.extend(score_data_to_insert)
                        ensembl_gene_id_to_canonical_mapping[gene_id] = True
                    elif canonical == "false":
                        ensembl_gene_id_to_canonical_mapping[gene_id] = False
                else:
                    if ensembl_gene_id_to_canonical_mapping[gene_id]:
                        if canonical == "true":
                            print(
                                f"Warning: Gene 'NA' with Gene ID '{gene_id}' has multiple canonical transcripts in input file data. Skipped gene."
                            )
                    else:
                        if canonical == "true":
                            gene = ensembl_to_gene_mapping[gene_id]
                            score_data_to_insert = validate_scores(
                                locus_to_locus_id_mapping,
                                gene,
                                gene_id,
                                pli_score,
                                loeuf_score,
                                pli_attrib_id,
                                loeuf_attrib_id,
                                source_id,
                                True,
                            )
                            data_to_insert.extend(score_data_to_insert)
                            ensembl_gene_id_to_canonical_mapping[gene_id] = True
            else:
                if gene_id not in not_found_ensembl_gene_id_list:
                    print(
                        f"Warning: Gene 'NA' with Gene ID '{gene_id}' not found in G2P DB. Skipped gene."
                    )
                    not_found_ensembl_gene_id_list.append(gene_id)
    for key, value in ensembl_gene_id_to_canonical_mapping.items():
        if not value:
            print(
                f"Warning: Gene 'NA' with Gene ID '{key}' does not have canonical transcript in input file data. Skipped gene."
            )
    return data_to_insert


def get_unique_gene_count(final_scores: list[tuple[str, str, str, str, str]]) -> int:
    """
    Calculate number of unique genes.

    Args:
        final_scores (list[tuple[str, str, str, str, str]]): list of tuples containing db values to insert

    Returns:
        int: number of unique genes
    """
    unique_gene_set = set()
    for item in final_scores:
        gene = item[0]
        unique_gene_set.add(gene)
    return len(unique_gene_set)


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
        dict[str, str]: dict of locus and associated ID
    """
    locus_to_locus_id_mapping = {}

    sql_get_locus = """ SELECT id, name from locus """
    sql_get_locus_attrib = """ SELECT locus_id, value from locus_attrib """

    database = MySQLdb.connect(
        host=db_host, port=db_port, user=user, passwd=password, db=db_name
    )
    cursor = database.cursor()

    cursor.execute(sql_get_locus)
    data = cursor.fetchall()
    for row in data:
        locus_id, locus = row
        if locus not in locus_to_locus_id_mapping:
            locus_to_locus_id_mapping[locus] = locus_id

    cursor.execute(sql_get_locus_attrib)
    data_2 = cursor.fetchall()
    for row in data_2:
        locus_id, locus = row
        if locus not in locus_to_locus_id_mapping:
            locus_to_locus_id_mapping[locus] = locus_id

    cursor.close()
    database.close()

    return locus_to_locus_id_mapping


def get_ensembl_id_from_g2p_db(
    db_host: str, db_port: int, db_name: str, user: str, password: str
) -> dict[str, str]:
    """
    Retrieves Ensembl IDs from the G2P database.

    Args:
        db_host (str): hostname
        db_port (int): port
        db_name (str): G2P database name
        user (str): username
        password (str): password

    Returns:
        dict[str, str]: dict of Ensembl ID and associated gene name
    """
    ensembl_to_gene_mapping = {}

    sql_get_ensembl_gene_id = """ SELECT identifier, name from locus_identifier join locus on locus_identifier.locus_id = locus.id """

    database = MySQLdb.connect(
        host=db_host, port=db_port, user=user, passwd=password, db=db_name
    )
    cursor = database.cursor()

    cursor.execute(sql_get_ensembl_gene_id)
    data = cursor.fetchall()
    for row in data:
        ensembl_id, gene = row
        if ensembl_id not in ensembl_to_gene_mapping:
            ensembl_to_gene_mapping[ensembl_id] = gene

    cursor.close()
    database.close()

    return ensembl_to_gene_mapping


def read_input_file(
    file: str,
    locus_to_locus_id_mapping: dict[str, str],
    ensembl_to_gene_mapping: dict[str, str],
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
        locus_to_locus_id_mapping(dict[str, str]): dict of genes and associated ids
        ensembl_to_gene_mapping(dict[str, str]): dict of ensembl gene id and associated gene names
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

    empty_genes_data_list = []
    with open(file, "r", encoding="utf-8") as f:
        previous_gene = ""
        previous_ensembl_gene_id = ""
        is_success = False
        # Below logic will process valid genes; 'NA' genes data will be processed afterwards
        for line in f:
            if not line.startswith("gene"):
                line_data = line.split("\t")
                current_gene = line_data[0]
                gene_id = line_data[1]
                canonical = line_data[3]
                pli_score = line_data[18]
                loeuf_score = line_data[22]
                if current_gene == "NA":
                    empty_genes_data_list.append(
                        (gene_id, canonical, pli_score, loeuf_score)
                    )
                else:
                    if current_gene != previous_gene:
                        validate_previous_gene_data(
                            locus_to_locus_id_mapping,
                            ensembl_to_gene_mapping,
                            previous_gene,
                            previous_ensembl_gene_id,
                            is_success,
                        )
                        previous_gene = current_gene
                        previous_ensembl_gene_id = ""
                        is_success = False
                    if gene_id.startswith("ENSG"):
                        previous_ensembl_gene_id = gene_id
                    if (
                        current_gene in locus_to_locus_id_mapping
                        and gene_id.startswith("ENSG")
                        and canonical == "true"
                    ):
                        if not is_success:
                            score_data_to_insert = validate_scores(
                                locus_to_locus_id_mapping,
                                current_gene,
                                gene_id,
                                pli_score,
                                loeuf_score,
                                pli_attrib_id,
                                loeuf_attrib_id,
                                source_id,
                            )
                            final_data_to_insert.extend(score_data_to_insert)
                            is_success = True
                        else:
                            print(
                                f"Warning: Gene '{current_gene}' has multiple canonical transcripts in input file data. Processed first match but skipped row with Gene '{current_gene}', Gene ID '{gene_id}' and canonical value '{canonical}'."
                            )
        validate_previous_gene_data(
            locus_to_locus_id_mapping,
            ensembl_to_gene_mapping,
            previous_gene,
            previous_ensembl_gene_id,
            is_success,
        )

    # Process 'NA' genes data
    empty_genes_data_to_insert = process_empty_genes(
        empty_genes_data_list,
        ensembl_to_gene_mapping,
        locus_to_locus_id_mapping,
        pli_attrib_id,
        loeuf_attrib_id,
        source_id,
    )
    final_data_to_insert.extend(empty_genes_data_to_insert)

    if len(final_data_to_insert) == 0:
        sys.exit("Error: No valid data found in input file.")

    return final_data_to_insert


def print_data_to_insert(final_scores: list[tuple[str, str, str, str, str]]) -> None:
    """
    Prints stats of the data to insert.

    Args:
        final_scores (list[tuple[str, str, str, str, str]]): list of tuples containing db values to insert

    Returns:
        None
    """
    print("----")
    print("Number of rows to be inserted: ", len(final_scores))
    print("Number of genes with scores: ", get_unique_gene_count(final_scores))
    print("----")


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
    get_attrib_id_query = """ SELECT id from attrib where value = %s """

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
    db_port = int(config["database"]["port"])
    db_name = config["database"]["name"]
    user = config["database"]["user"]
    pwd = config["database"]["password"]

    print("Validating input file...")
    validate_input_file(file)
    print("Validating input file... done\n")

    print("Getting locus id from G2P DB...")
    locus_to_locus_id_mapping = get_locus_id_from_g2p_db(
        db_host, db_port, db_name, user, pwd
    )
    print("Getting locus id from G2P DB... done\n")

    print("Getting ensembl gene id from G2P DB...")
    ensembl_to_gene_mapping = get_ensembl_id_from_g2p_db(
        db_host, db_port, db_name, user, pwd
    )
    print("Getting ensembl gene id from G2P DB... done\n")

    print("Reading input file...")
    final_scores = read_input_file(
        file,
        locus_to_locus_id_mapping,
        ensembl_to_gene_mapping,
        db_host,
        db_port,
        db_name,
        user,
        pwd,
    )
    print("Reading input file... done\n")

    print_data_to_insert(final_scores)

    print("Inserting data into gene_stats table...")
    insert_into_gene_stats(final_scores, db_host, db_port, db_name, user, pwd)
    print("Inserting data into gene_stats table... done\n")

    print("Inserting data into meta table...")
    insert_details_into_meta(db_host, db_port, db_name, user, pwd)
    print("Inserting data into meta table... done\n")


if __name__ == "__main__":
    main()
