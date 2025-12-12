#!/usr/bin/env python3

import argparse
import configparser
import gzip
import os
import re
import sys
import urllib.request
from datetime import datetime

import MySQLdb


"""
    Script to update genes in G2P database based on a new Ensembl GTF file.
    By default, it updates/adds genes based on the Ensembl stable IDs.
    To only update the gene symbols based on HGNC data, use the --only_update_gene_symbols option.

    Options:
            --config:       Config file with details to the G2P database (mandatory)
                            File format is the following:
                                [g2p_database]
                                host = <>
                                port = <>
                                user = <>
                                password = <>
                                name = <>
            --working_dir:  Working directory to store report files (mandatory)
            --version:      Ensembl version only required if going to run the genes update (optional)
            --only_update_gene_symbols:  Only update gene symbols based on HGNC data (optional)
"""


def locus_id_foreign_key_check(
    db_host: str,
    db_port: int,
    db_name: str,
    user: str,
    password: str,
    g2p_tables_with_locus_id_link: list,
) -> None:
    """
    Run a FK check on the tables that have a link to the locus table.

    Args:
        db_host (str): G2P host name
        db_port (int): port number
        db_name (str): G2P database name
        user (str): user
        password (str): password
        g2p_tables_with_locus_id_link (list): list of G2P tables with a link to the locus table
    """

    db = MySQLdb.connect(
        host=db_host, port=db_port, user=user, passwd=password, db=db_name
    )
    cursor = db.cursor()

    for table in g2p_tables_with_locus_id_link:
        column_name = "gene_id"

        if table.startswith("locus_") or table.startswith("gene2phenotype_app_"):
            column_name = "locus_id"

        sql = f"""
                    SELECT COUNT(t.id) FROM {table} t
                    LEFT JOIN locus l ON l.id = t.{column_name}
                    WHERE l.id IS NULL AND t.{column_name} IS NOT NULL
                """

        cursor.execute(sql)
        data = cursor.fetchone()[0]
        if data != 0:
            sys.exit(
                f"Found {data} unlinked entries in table {table} after locus_id foreign key check\n"
            )

    db.close()


def check_disease_names(
    working_dir: str, db_host: str, db_port: int, db_name: str, user: str, password: str
) -> None:
    """
    Check which disease names do not include the latest gene symbol in the name.
    It writes the list of diseases to a file.

    Args:
        working_dir (str): working directory
        db_host (str): G2P host name
        db_port (int): port number
        db_name (str): G2P database name
        user (str): user
        password (str): password
    """
    report_file = os.path.join(working_dir, "report_diseases_to_update.txt")

    sql = """   SELECT g.stable_id,l.name,d.name from locus_genotype_disease lgd
                LEFT JOIN disease d ON d.id = lgd.disease_id
                LEFT JOIN g2p_stableid g ON g.id = lgd.stable_id
                LEFT JOIN locus l ON l.id = lgd.locus_id
                WHERE g.is_deleted = 0
          """

    db = MySQLdb.connect(
        host=db_host, port=db_port, user=user, passwd=password, db=db_name
    )
    cursor = db.cursor()
    cursor.execute(sql)
    data = cursor.fetchall()

    with open(report_file, "w") as wr:
        for row in data:
            if not (row[2].startswith(f"{row[1]}-") or row[2].startswith(f"{row[1]} ")):
                wr.write(
                    f"Disease to update: {row[0]}; gene '{row[1]}'; disease '{row[2]}'\n"
                )

    db.close()


def read_from_gtf(
    working_dir: str, ensembl_gtf: str, exclude_biotypes: list
) -> dict[str, dict]:
    """
    Method to read the Ensembl GTF file and create a list if genes to import/update.

    Args:
        working_dir (str): working directory
        ensembl_gtf (str): GTF file from Ensembl
        exclude_biotypes (list): list of biotypes to exclude from the update

    Returns:
        dict[str, dict]: dictionary of genes to import/update
    """
    gene_symbol_2_stable_id = {}
    gene_symbol_details = {}
    unique_gene_symbol_2_stable_id = {}
    unique_stable_id_2_gene_symbol = {}

    with gzip.open(ensembl_gtf, "rt", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue

            data = line.strip().split("\t")
            attribs_list = {}
            for attrib in data[8].split(";"):
                if (
                    attrib.strip().startswith("gene_id")
                    or attrib.strip().startswith("gene_name")
                    or attrib.strip().startswith("gene_biotype")
                ):
                    get_attrib_values = attrib.strip().replace('"', "").split(" ")
                    if get_attrib_values and len(get_attrib_values) == 2:
                        attribs_list[get_attrib_values[0]] = get_attrib_values[1]

            if (
                "gene_id" in attribs_list
                and attribs_list["gene_id"]
                and "gene_name" in attribs_list
                and attribs_list["gene_name"]
                and "gene_biotype" in attribs_list
                and attribs_list["gene_biotype"]
                and not any(
                    re.search(sub, attribs_list["gene_biotype"])
                    for sub in exclude_biotypes
                )
            ):
                try:
                    gene_symbol_2_stable_id[attribs_list["gene_name"]]
                except KeyError:
                    gene_symbol_2_stable_id[attribs_list["gene_name"]] = set()
                    gene_symbol_2_stable_id[attribs_list["gene_name"]].add(
                        attribs_list["gene_id"]
                    )
                    # Save the chr details for each gene symbol
                    gene_symbol_details[attribs_list["gene_name"]] = {
                        "chr": data[0],
                        "start": data[3],
                        "end": data[4],
                        "strand": data[6],
                    }
                else:
                    start = int(data[3])
                    end = int(data[4])
                    # Genes in YPAR regions have the same gene symbol but different gene ID
                    # These are valid genes
                    if not (
                        (
                            data[0] == "X"
                            and (
                                (start >= 10001 and end <= 2781479)
                                or (start >= 155701383 and end <= 156030895)
                            )
                        )
                        or (
                            data[0] == "Y"
                            and (
                                (start >= 10001 and end <= 2781479)
                                or (start >= 56887903 and end <= 57217415)
                            )
                        )
                    ):
                        gene_symbol_2_stable_id[attribs_list["gene_name"]].add(
                            attribs_list["gene_id"]
                        )

    # Write report files
    genes_output_file = os.path.join(working_dir, "ensembl_genes_grch38.txt")
    error_log_file = os.path.join(working_dir, "ensembl_genes_grch38_error.log")

    with open(error_log_file, "w") as wr, open(genes_output_file, "w") as wr_genes:
        for gene_symbol in gene_symbol_2_stable_id:
            if len(gene_symbol_2_stable_id[gene_symbol]) > 1:
                wr.write(
                    f"ERROR: more than one stable_id for {gene_symbol}: {gene_symbol_2_stable_id[gene_symbol]}\n"
                )
            else:
                unique_gene_symbol_2_stable_id[gene_symbol] = list(
                    gene_symbol_2_stable_id[gene_symbol]
                )[0]
                wr_genes.write(
                    f"{unique_gene_symbol_2_stable_id[gene_symbol]}\t{gene_symbol}\n"
                )
                # save the ensembl id as key for easier comparison later
                unique_stable_id_2_gene_symbol[
                    unique_gene_symbol_2_stable_id[gene_symbol]
                ] = {
                    "gene_symbol": gene_symbol,
                    "chr": gene_symbol_details[gene_symbol]["chr"],
                    "start": gene_symbol_details[gene_symbol]["start"],
                    "end": gene_symbol_details[gene_symbol]["end"],
                    "strand": gene_symbol_details[gene_symbol]["strand"],
                }

    return unique_stable_id_2_gene_symbol


def get_g2p_genes(
    db_host: str, db_port: int, db_name: str, user: str, password: str
) -> tuple[dict[str, dict], dict[str, str]]:
    """
    Method to fetch all the genes stored in the G2P database.
    It returns the same information in two separate dictionaries with different formats:
        g2p_genes: complete dictionary of genes data where the key is the symbol
        g2p_genes_by_symbol: dictionary of genes where the key is the symbol and value is the Ensembl ID

    Args:
        db_host (str): G2P host name
        db_port (int): port number
        db_name (str): G2P database name
        user (str): user
        password (str): password

    Returns:
        tuple[dict[str, dict], dict[str, str]]: two dictionaries with the G2P genes data
    """
    g2p_genes = {}
    g2p_genes_by_symbol = {}

    sql = """
                SELECT l.name, li.identifier, l.id, l.start, l.end, seq.name FROM locus l
                LEFT JOIN locus_identifier li ON li.locus_id = l.id
                LEFT JOIN source s ON s.id = li.source_id
                LEFT JOIN sequence seq ON seq.id = l.sequence_id
                WHERE s.name = 'Ensembl'
            """

    db = MySQLdb.connect(
        host=db_host, port=db_port, user=user, passwd=password, db=db_name
    )
    cursor = db.cursor()
    cursor.execute(sql)
    data = cursor.fetchall()
    for row in data:
        g2p_genes[row[1]] = {
            "gene_symbol": row[0],
            "locus_id": row[2],
            "locus_start": row[3],
            "locus_end": row[4],
            "locus_chr": row[5],
        }
        g2p_genes_by_symbol[row[0]] = row[1]

    db.close()

    return g2p_genes, g2p_genes_by_symbol


def get_g2p_genes_hgnc(
    db_host: str, db_port: int, db_name: str, user: str, password: str
) -> tuple[dict[str, dict], dict[str, dict]]:
    """
    Method to fetch the full G2P gene set with the HGNC and OMIM IDs.
    Returns a dictionary where the key is the Ensembl gene ID.

    Args:
        db_host (str): G2P host name
        db_port (int): port number
        db_name (str): G2P database name
        user (str): user
        password (str): password

    Returns two dictionaries:
        - g2p_genes_by_ensembl_id: keyed by Ensembl gene ID, with gene details and synonyms.
        - g2p_genes_by_symbol: keyed by gene symbol, with HGNC and OMIM IDs.

    Called by: update_xrefs()
    """
    g2p_genes_by_ensembl_id = {}  # key = Ensembl gene ID
    g2p_genes_by_symbol = {}  # key = gene symbol

    sql = """
                SELECT l.name, li.identifier, l.id, la.value FROM locus l
                LEFT JOIN locus_identifier li ON li.locus_id = l.id
                LEFT JOIN locus_attrib la ON la.locus_id = l.id
                LEFT JOIN source s ON s.id = li.source_id
                WHERE s.name = 'Ensembl'
            """

    sql_ids = """
                    SELECT l.name, li.identifier FROM locus l
                    LEFT JOIN locus_identifier li ON li.locus_id = l.id
                """

    db = MySQLdb.connect(
        host=db_host, port=db_port, user=user, passwd=password, db=db_name
    )
    cursor = db.cursor()
    cursor.execute(sql)
    data = cursor.fetchall()
    for row in data:
        try:
            g2p_genes_by_ensembl_id[row[1]]
        except KeyError:
            g2p_genes_by_ensembl_id[row[1]] = {
                "gene_symbol": row[0],
                "locus_id": row[2],
                "synonyms": [row[3]],
            }
        else:
            g2p_genes_by_ensembl_id[row[1]]["synonyms"].append(row[3])

    cursor.execute(sql_ids)
    data_ids = cursor.fetchall()
    for row in data_ids:
        if row[0] not in g2p_genes_by_symbol:
            g2p_genes_by_symbol[row[0]] = {}
            if row[1].startswith("HGNC"):
                g2p_genes_by_symbol[row[0]]["hgnc_id"] = row[1]
            elif row[1].isdigit():
                g2p_genes_by_symbol[row[0]]["omim_id"] = row[1]
        else:
            if row[1].startswith("HGNC"):
                g2p_genes_by_symbol[row[0]]["hgnc_id"] = row[1]
            elif row[1].isdigit():
                g2p_genes_by_symbol[row[0]]["omim_id"] = row[1]

    db.close()

    return g2p_genes_by_ensembl_id, g2p_genes_by_symbol


def get_g2p_genes_not_used(
    db_host: str, db_port: int, db_name: str, user: str, password: str
) -> dict[str, str]:
    """
    Method to fetch a list of G2P genes that are not linked to any tables.

    Args:
        db_host (str): G2P host name
        db_port (int): port number
        db_name (str): G2P database name
        user (str): user
        password (str): password

    Returns:
        dict[str, str]: List of G2P genes not used in other tables
    """

    g2p_genes_not_used = {}

    sql_not_used = """
                        SELECT l.name, li.identifier, l.id FROM locus l
                        LEFT JOIN locus_identifier li ON l.id = li.locus_id
                        LEFT JOIN source s ON s.id = li.source_id
                        LEFT JOIN locus_genotype_disease lgd ON lgd.locus_id = l.id
                        LEFT JOIN gene_stats g ON g.gene_id = l.id
                        LEFT JOIN uniprot_annotation u ON u.gene_id = l.id
                        LEFT JOIN gene_disease d ON d.gene_id = l.id
                        WHERE l.id is not null and lgd.locus_id is null and g.gene_id is null and 
                        u.gene_id is null and d.gene_id is null and s.name = 'Ensembl'
                    """

    db = MySQLdb.connect(
        host=db_host, port=db_port, user=user, passwd=password, db=db_name
    )
    cursor = db.cursor()
    cursor.execute(sql_not_used)
    data_other = cursor.fetchall()
    for row in data_other:
        g2p_genes_not_used[row[1]] = {"name": row[0], "id": row[2]}

    db.close()

    return g2p_genes_not_used


def update_genes(
    working_dir: str,
    g2p_gene_ids: dict[str, dict],
    g2p_genes_by_symbol: dict[str, str],
    unique_stable_id_2_gene_symbol: dict[str, dict],
    db_host: str,
    db_port: int,
    db_name: str,
    user: str,
    password: str,
) -> None:
    """
    Method to run the genes update.
    It writes two report files:
        report_gene_updates.txt
        report_new_genes.txt

    Args:
        working_dir (str): working directory
        g2p_gene_ids (dict[str, dict]): complete dictionary of genes data where the key is the symbol
        g2p_genes_by_symbol (dict[str, str]): dictionary of genes where the key is the symbol and value is the Ensembl ID
        unique_stable_id_2_gene_symbol (dict[str, dict]): dictionary of genes to import/update
        db_host (str): G2P host name
        db_port (int): port number
        db_name (str): G2P database name
        user (str): user
        password (str): password
    """

    update_report = os.path.join(working_dir, "report_gene_updates.txt")
    new_genes_report = os.path.join(working_dir, "report_new_genes.txt")

    sequence_chrs = {}

    sql_sel_sequences = """
                            SELECT id, name FROM sequence
                        """

    sql_sel_attrib_gene = """
                                SELECT id FROM attrib WHERE value = 'gene'
                            """

    sql_sel_source = """
                            SELECT id FROM source WHERE name = 'Ensembl'
                        """

    sql_insert_identifier = """
                                INSERT INTO locus_identifier(identifier, description, locus_id, source_id)
                                VALUES(%s, %s, %s, %s)
                            """

    sql_update_stable_id = """
                                UPDATE locus_identifier SET identifier = %s
                                WHERE source_id = %s AND locus_id = %s
                            """

    sql_insert_gene = """
                            INSERT INTO locus(start, end, strand, name, type_id, sequence_id)
                            VALUES(%s, %s, %s, %s, %s, %s)
                        """

    db = MySQLdb.connect(
        host=db_host, port=db_port, user=user, passwd=password, db=db_name
    )
    cursor = db.cursor()
    # Get all sequences (chrs)
    cursor.execute(sql_sel_sequences)
    data = cursor.fetchall()
    for row in data:
        sequence_chrs[row[1]] = row[0]

    cursor.execute(sql_sel_attrib_gene)
    gene_attrib_id = cursor.fetchone()
    cursor.execute(sql_sel_source)
    source_id = cursor.fetchone()

    with open(update_report, "w") as wr, open(new_genes_report, "w") as wr_new:
        for stable_id in unique_stable_id_2_gene_symbol:
            new_gene_symbol = unique_stable_id_2_gene_symbol[stable_id]["gene_symbol"]

            try:
                g2p_gene_ids[stable_id]["gene_symbol"]
            except KeyError:
                # The stable_id is not found in G2P, it could mean one of the following:
                # 1) gene is in g2p but in the gtf file the stable id has been updated for the gene symbol
                # 2) this is a new gene
                if new_gene_symbol in g2p_genes_by_symbol:
                    # Scenario 1
                    g2p_current_stable_id = g2p_genes_by_symbol[new_gene_symbol]
                    cursor.execute(
                        sql_update_stable_id,
                        [
                            stable_id,
                            source_id,
                            g2p_gene_ids[g2p_current_stable_id]["locus_id"],
                        ],
                    )
                    db.commit()
                    # Write to report
                    wr.write(
                        f"UPDATE: locus_id = {g2p_gene_ids[g2p_current_stable_id]['locus_id']} gene symbol {new_gene_symbol} new stable_id {stable_id} (previous stable_id {g2p_genes_by_symbol[new_gene_symbol]})\n"
                    )
                else:
                    # Scenario 2
                    if unique_stable_id_2_gene_symbol[stable_id]["strand"] == "-":
                        strand = -1
                    else:
                        strand = 1

                    # Add gene
                    cursor.execute(
                        sql_insert_gene,
                        [
                            unique_stable_id_2_gene_symbol[stable_id]["start"],
                            unique_stable_id_2_gene_symbol[stable_id]["end"],
                            strand,
                            new_gene_symbol,
                            gene_attrib_id,
                            sequence_chrs[
                                unique_stable_id_2_gene_symbol[stable_id]["chr"]
                            ],
                        ],
                    )
                    db.commit()
                    locus_id = cursor.lastrowid
                    # Add Ensembl stable_id
                    cursor.execute(
                        sql_insert_identifier, [stable_id, None, locus_id, source_id]
                    )
                    db.commit()
                    # Write to report
                    wr_new.write(
                        f"ADD: new gene {new_gene_symbol} stable_id {stable_id}\n"
                    )

    db.close()


def looks_like_identifier(symbol: str) -> bool:
    """
    Method to check if a string is a gene identifier.
    """
    return bool(re.match(r"^[A-Z]+[0-9]+\.[0-9]+", symbol))


def delete_outdated_locus(
    working_dir: str,
    db_host: str,
    db_port: int,
    db_name: str,
    user: str,
    password: str,
    g2p_gene_ids: dict[str, dict],
    unique_stable_id_2_gene_symbol: dict[str, dict],
) -> None:
    """
    Method to delete the genes that are outdated, e.g. not found in the new GTF file.

    Args:
        working_dir (str): working directory
        db_host (str): G2P host name
        db_port (int): port number
        db_name (str): G2P database name
        user (str): user
        password (str): password
        g2p_gene_ids (dict[str, dict]): complete dictionary of genes data where the key is the symbol
        unique_stable_id_2_gene_symbol (dict[str, dict]): dictionary of genes to import/update
    """
    # Prepare output file
    report = os.path.join(working_dir, "report_outdated_genes.txt")

    sql_locus = """
                    DELETE FROM locus WHERE id = %s
                """

    sql_locus_1 = """
                        DELETE FROM locus_attrib WHERE locus_id = %s
                    """

    sql_locus_2 = """
                        DELETE FROM locus_identifier WHERE locus_id = %s
                    """

    db = MySQLdb.connect(
        host=db_host, port=db_port, user=user, passwd=password, db=db_name
    )
    cursor = db.cursor()

    # Fetch the list of genes that are not being used in any table
    # The dictionary key is the gene stable id
    g2p_genes_not_used = get_g2p_genes_not_used(
        db_host, db_port, db_name, user, password
    )

    with open(report, "w") as wr:
        for stable_id in g2p_gene_ids:
            if stable_id not in unique_stable_id_2_gene_symbol:
                gene_symbol = g2p_gene_ids[stable_id]["gene_symbol"]
                if stable_id in g2p_genes_not_used:
                    locus_id = g2p_genes_not_used[stable_id]["id"]
                    cursor.execute(sql_locus_2, [locus_id])
                    cursor.execute(sql_locus_1, [locus_id])
                    cursor.execute(sql_locus, [locus_id])
                    db.commit()
                    wr.write(
                        f"INFO: outdated {stable_id} ({gene_symbol}) deleted from G2P\n"
                    )
                else:
                    wr.write(
                        f"WARNING: outdated locus used in G2P {stable_id} ({gene_symbol})\n"
                    )

    cursor.close()
    db.close()


def update_xrefs(
    working_dir: str,
    hgnc_file: str,
    db_host: str,
    db_port: int,
    db_name: str,
    user: str,
    password: str,
) -> None:
    """
    Main method to update the genes identifiers and gene symbols.
    This method calls:
        get_g2p_genes_hgnc()
        update_g2p_identifier()
        update_locus_name()
    It writes a report to report_hgnc_updates.txt

    Args:
        working_dir (str): working directory
        hgnc_file (str): File containing the HGNC data
        db_host (str): G2P host name
        db_port (int): port number
        db_name (str): G2P database name
        user (str): user
        password (str): password
    """
    report_hgnc_file = os.path.join(working_dir, "report_hgnc_updates.txt")

    # Get a new dump of the g2p genes - after they were updated with the gft file
    g2p_genes_by_ensembl_id, g2p_genes_by_symbol = get_g2p_genes_hgnc(
        db_host, db_port, db_name, user, password
    )

    mappings = {}

    with open(
        hgnc_file,
        "r",
    ) as fh:
        for line in fh:
            if line.startswith("hgnc_id"):
                header = line.strip().split("\t")
                continue
            elif not line.startswith("HGNC:"):
                continue

            data = line.strip().split("\t")
            data = dict(zip(header, data))
            hgnc_gene_symbol = data["symbol"]

            if "ensembl_gene_id" not in data:
                continue
            else:
                ensembl_gene_id = data["ensembl_gene_id"]

            for id_type in ["hgnc_id", "prev_symbol", "omim_id"]:
                try:
                    value = data[id_type]
                except KeyError:
                    continue
                else:
                    # clean the value before saving it
                    value = value.replace('"', "")
                    if ensembl_gene_id not in mappings:
                        mappings[ensembl_gene_id] = {}
                        mappings[ensembl_gene_id][id_type] = [value]
                    elif id_type in mappings[ensembl_gene_id]:
                        mappings[ensembl_gene_id][id_type].append(value)
                    else:
                        mappings[ensembl_gene_id][id_type] = [value]
            # Plus gene symbol
            if ensembl_gene_id in mappings:
                mappings[ensembl_gene_id]["gene_symbol"] = hgnc_gene_symbol

    # Update
    with open(
        report_hgnc_file,
        "w",
    ) as wr:
        for ensembl_gene_id in mappings:
            # Get current g2p data for this gene ensembl id
            try:
                g2p_data = g2p_genes_by_ensembl_id[ensembl_gene_id]
            except KeyError:
                continue
            else:
                current_gene_symbol = g2p_data["gene_symbol"]

                # Get the identifiers from G2P for this gene
                current_identifiers = g2p_genes_by_symbol[current_gene_symbol]

                # Update or add HGNC ID
                if len(mappings[ensembl_gene_id]["hgnc_id"]) == 1:
                    if (
                        "hgnc_id" in current_identifiers
                        and mappings[ensembl_gene_id]["hgnc_id"][0]
                        != current_identifiers["hgnc_id"]
                    ):
                        update_g2p_identifier(
                            "HGNC_update",
                            mappings[ensembl_gene_id]["hgnc_id"][0],
                            g2p_data["locus_id"],
                            db_host,
                            db_port,
                            db_name,
                            user,
                            password,
                        )
                        wr.write(
                            f"UPDATE HGNC ID: locus_id = {g2p_data['locus_id']} gene symbol {current_gene_symbol} {mappings[ensembl_gene_id]['hgnc_id'][0]}\n"
                        )
                    elif "hgnc_id" not in current_identifiers:
                        update_g2p_identifier(
                            "HGNC_insert",
                            mappings[ensembl_gene_id]["hgnc_id"][0],
                            g2p_data["locus_id"],
                            db_host,
                            db_port,
                            db_name,
                            user,
                            password,
                        )
                        wr.write(
                            f"ADD HGNC ID: locus_id = {g2p_data['locus_id']} gene symbol {current_gene_symbol} {mappings[ensembl_gene_id]['hgnc_id'][0]}\n"
                        )

                # Update or add OMIM ID
                if "omim_id" in mappings[ensembl_gene_id]:
                    # Format OMIM IDs
                    omim_ids = mappings[ensembl_gene_id]["omim_id"][0].split("|")
                    if (
                        len(omim_ids) == 1
                        and omim_ids[0] != ""
                        and "omim_id" in current_identifiers
                        and omim_ids[0] != current_identifiers["omim_id"]
                    ):
                        update_g2p_identifier(
                            "OMIM_update",
                            omim_ids[0],
                            g2p_data["locus_id"],
                            db_host,
                            db_port,
                            db_name,
                            user,
                            password,
                        )
                        wr.write(
                            f"UPDATE OMIM ID: locus_id = {g2p_data['locus_id']} gene symbol {current_gene_symbol} {omim_ids[0]}\n"
                        )
                    elif (
                        "omim_id" not in current_identifiers
                        and omim_ids[0] != ""
                        and len(omim_ids) == 1
                    ):
                        update_g2p_identifier(
                            "OMIM_insert",
                            omim_ids[0],
                            g2p_data["locus_id"],
                            db_host,
                            db_port,
                            db_name,
                            user,
                            password,
                        )
                        wr.write(
                            f"ADD OMIM ID: locus_id = {g2p_data['locus_id']} gene symbol {current_gene_symbol} {omim_ids[0]}\n"
                        )

                # Update gene symbol
                new_gene_symbol = mappings[ensembl_gene_id]["gene_symbol"]
                add_synonym = False
                if (
                    "synonyms" in g2p_data
                    and current_gene_symbol not in g2p_data["synonyms"]
                ):
                    add_synonym = True
                if current_gene_symbol != new_gene_symbol:
                    update_locus_name(
                        new_gene_symbol,
                        current_gene_symbol,
                        g2p_data["locus_id"],
                        add_synonym,
                        db_host,
                        db_port,
                        db_name,
                        user,
                        password,
                    )
                    wr.write(
                        f"UPDATE GENE SYMBOL: locus_id = {g2p_data['locus_id']} old gene symbol '{current_gene_symbol}' new gene symbol '{new_gene_symbol}'\n"
                    )


def update_g2p_identifier(
    type: str,
    new_hgnc_id: str,
    locus_id: int,
    db_host: str,
    db_port: int,
    db_name: str,
    user: str,
    password: str,
) -> None:
    """
    Method to run the sql queries to update the gene identifiers.

    Args:
        type (str): HGNC or OMIM
        new_hgnc_id (str): new HGNC ID
        locus_id (int): locus ID that is going to be updated
        db_host (str): G2P host name
        db_port (int): port number
        db_name (str): G2P database name
        user (str): user
        password (str): password
    """

    source_name, query_type = type.split("_")

    sql_source = """
                        SELECT id FROM source
                        WHERE name = %s
                    """

    sql_update_hgnc = """
                            UPDATE locus_identifier SET identifier = %s
                            WHERE locus_id = %s AND source_id = %s
                        """

    sql_insert_hgnc = """
                            INSERT INTO locus_identifier(identifier, description, locus_id, source_id)
                            VALUES(%s, %s, %s, %s)
                        """

    db = MySQLdb.connect(
        host=db_host, port=db_port, user=user, passwd=password, db=db_name
    )
    cursor = db.cursor()
    cursor.execute(sql_source, [source_name])
    source_id = cursor.fetchone()

    if not source_id:
        print(f"ERROR: cannot find source name {source_name} in g2p")

    if query_type == "update":
        cursor.execute(sql_update_hgnc, [new_hgnc_id, locus_id, source_id])
    else:
        cursor.execute(sql_insert_hgnc, [new_hgnc_id, None, locus_id, source_id])
    db.commit()
    db.close()


def update_locus_name(
    symbol: str,
    old_symbol: str,
    locus_id: int,
    add_synonym: bool,
    db_host: str,
    db_port: int,
    db_name: str,
    user: str,
    password: str,
) -> None:
    """
    Method to run the sql query to update the gene symbols in table 'locus'.

    Args:
        symbol (str): new gene symbol
        old_symbol (str): old gene symbol
        locus_id (int): locus ID that is going to be updated
        add_synonym (bool): flag if synonym should be added
        db_host (str): G2P host name
        db_port (int): port number
        db_name (str): G2P database name
        user (str): user
        password (str): password
    """
    sql_sel_attribs = "SELECT id FROM attrib_type WHERE code = 'gene_synonym'"

    sql_sel_source = "SELECT id FROM source WHERE name = 'Ensembl'"

    sql_update_gene = """UPDATE locus SET name = %s WHERE id = %s"""

    sql_insert_syn = """
                        INSERT INTO locus_attrib(value, is_deleted, attrib_type_id, locus_id, source_id)
                        VALUES(%s, %s, %s, %s, %s)
                    """

    db = MySQLdb.connect(
        host=db_host, port=db_port, user=user, passwd=password, db=db_name
    )
    cursor = db.cursor()
    cursor.execute(sql_sel_attribs)
    attrib_type_id = cursor.fetchone()
    cursor.execute(sql_sel_source)
    source_id = cursor.fetchone()

    # Update gene symbol
    cursor.execute(sql_update_gene, [symbol, locus_id])

    # Add old gene symbol as a synonym (locus_attrib)
    if add_synonym:
        cursor.execute(
            sql_insert_syn, [old_symbol, 0, attrib_type_id, locus_id, source_id]
        )

    db.commit()
    db.close()


def update_meta(
    db_host: str,
    db_port: int,
    db_name: str,
    user: str,
    password: str,
    version: str,
    update_gene_symbol: bool,
) -> None:
    """
    Method to add info about update to the meta table.

    Args:
        db_host (str): G2P host name
        db_port (int): port number
        db_name (str): G2P database name
        user (str): user
        password (str): password
        version (str): Ensembl version
        update_gene_symbol (bool): update only the gene symbols
    """
    date_now = datetime.now()
    type_of_update = "locus_gene_update"

    # The script only updated the gene symbols
    if update_gene_symbol:
        type_of_update = "locus_gene_symbol_update"
        version = date_now.strftime("%Y-%m-%d")
        description = "Update gene symbols from HGNC"
        sql_version = "SELECT id FROM source WHERE name = 'HGNC'"
    else:
        # The script updated the whole gene set (including gene symbols)
        description = f"Update genes to Ensembl release {version}"
        sql_version = "SELECT id FROM source WHERE name = 'Ensembl'"

    sql = """
                INSERT INTO meta(`key`, date_update, is_public, description, version, source_id)
                VALUES(%s, %s, %s, %s, %s, %s)
            """

    db = MySQLdb.connect(
        host=db_host, port=db_port, user=user, passwd=password, db=db_name
    )
    cursor = db.cursor()
    # Fetch the source_id
    cursor.execute(sql_version)
    source_id = cursor.fetchone()
    # Add new row to meta
    cursor.execute(
        sql,
        [type_of_update, date_now, 0, description, version, source_id[0]],
    )
    db.commit()
    db.close()


def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--config", required=True, help="Config file")
    parser.add_argument(
        "--version", required=False, default=None, help="Ensembl release version"
    )
    parser.add_argument("--working_dir", required=True, help="Working directory")
    parser.add_argument(
        "--only_update_gene_symbol",
        action="store_true",
        help="Only update the gene symbol and identifiers from the HGNC file",
    )
    args = parser.parse_args()

    config_file = args.config
    version = args.version
    working_dir = args.working_dir
    only_update_gene_symbol = args.only_update_gene_symbol

    if not only_update_gene_symbol and not version:
        sys.exit(
            "You have to provide an Ensembl version (--version) or run with option --only_update_gene_symbol"
        )

    hgnc_file_url = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
    hgnc_file = os.path.join(working_dir, "hgnc_complete_set.txt")

    if not os.path.exists(working_dir):
        sys.exit(f"Invalid directory '{working_dir}'")

    try:
        urllib.request.urlretrieve(hgnc_file_url, hgnc_file)
    except (urllib.error.URLError, urllib.error.HTTPError) as error:
        sys.exit(f"Problem while fetching HGNC file '{hgnc_file_url}': {error}")

    # Load the config file
    config = configparser.ConfigParser()
    config.read(config_file)

    try:
        g2p_db_details = config["g2p_database"]
    except KeyError:
        sys.exit("Config: [g2p_database] is missing from the config file")
    else:
        db_host = g2p_db_details["host"]
        db_port = int(g2p_db_details["port"])
        db_name = g2p_db_details["name"]
        user = g2p_db_details["user"]
        password = g2p_db_details["password"]

    exclude_biotypes = ["pseudogene", "misc_RNA"]

    g2p_tables_with_locus_id_link = [
        "locus_identifier",
        "locus_attrib",
        "locus_genotype_disease",
        "uniprot_annotation",
        "gene_stats",
        "gene_disease",
        "gene2phenotype_app_historicallocusgenotypedisease",  # history table
    ]

    if not only_update_gene_symbol:
        print("Updating gene set...")
        # Download Ensembl file
        ensembl_gtf_url = f"https://ftp.ensembl.org/pub/release-{version}/gtf/homo_sapiens/Homo_sapiens.GRCh38.{version}.chr.gtf.gz"
        ensembl_gtf = os.path.join(
            working_dir, f"Homo_sapiens.GRCh38.{version}.chr.gtf.gz"
        )

        try:
            urllib.request.urlretrieve(ensembl_gtf_url, ensembl_gtf)
        except (urllib.error.URLError, urllib.error.HTTPError) as error:
            sys.exit(
                f"Problem while fetching Ensembl GTF file '{ensembl_gtf_url}': {error}"
            )

        # Check each table with locus id if all data in that column exists in the locus table
        locus_id_foreign_key_check(
            db_host, db_port, db_name, user, password, g2p_tables_with_locus_id_link
        )

        # Get the current G2P genes (Ensembl ID and corresponding gene symbol) -> g2p_gene_ids
        # Get the current G2P genes (gene symbol and corresponding Ensembl ID) -> g2p_genes_by_symbol
        g2p_gene_ids, g2p_genes_by_symbol = get_g2p_genes(
            db_host, db_port, db_name, user, password
        )

        # Get the Ensembl genes from the gtf file
        # Only considers genes that are linked to one Ensembl ID expect genes overlaping a YPAR region
        unique_stable_id_2_gene_symbol = read_from_gtf(
            working_dir, ensembl_gtf, exclude_biotypes
        )

        # If gene not in G2P insert new gene
        # This method does not update gene_symbol - this update is done later in method update_xrefs()
        # If the gene_symbol is already present in the locus table only update the ensembl stable_id
        update_genes(
            working_dir,
            g2p_gene_ids,
            g2p_genes_by_symbol,
            unique_stable_id_2_gene_symbol,
            db_host,
            db_port,
            db_name,
            user,
            password,
        )

        # Check if there are G2P genes that should be removed because they are not valid in the gtf file
        delete_outdated_locus(
            working_dir,
            db_host,
            db_port,
            db_name,
            user,
            password,
            g2p_gene_ids,
            unique_stable_id_2_gene_symbol,
        )
        print("Updating gene set... done")

    # Update the HGNC IDs and the gene symbols
    # This update uses the Ensembl gene ID to compare the genes
    # We need a new dump of the genes (after all the previous updates)
    print("Updating gene symbols and identifiers...")
    update_xrefs(working_dir, hgnc_file, db_host, db_port, db_name, user, password)
    print("Updating gene symbols and identifiers... done")

    print("Running checks...")
    # Check which disease names have to be updated
    check_disease_names(working_dir, db_host, db_port, db_name, user, password)

    # Check again each table with locus id if all data in that column exists in the locus table
    locus_id_foreign_key_check(
        db_host, db_port, db_name, user, password, g2p_tables_with_locus_id_link
    )
    print("Running checks... done")

    # Update meta
    update_meta(
        db_host, db_port, db_name, user, password, version, only_update_gene_symbol
    )


if __name__ == "__main__":
    main()
