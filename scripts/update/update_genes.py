#!/usr/bin/env python3

import os.path
import sys
import argparse
import MySQLdb
import requests
import urllib.request
import configparser
import gzip
import re

def locus_id_foreign_key_check(db_host, db_port, db_name, user, password, g2p_tables_with_locus_id_link):
    db = MySQLdb.connect(host=db_host, port=db_port, user=user, passwd=password, db=db_name)
    cursor = db.cursor()

    for table in g2p_tables_with_locus_id_link:
        column_name = "gene_id"

        if table.startswith("locus_") or table.startswith("gene2phenotype_app_"):
            column_name = "locus_id"

        sql =   f"""
                    SELECT COUNT(t.id) FROM {table} t
                    LEFT JOIN locus l ON l.id = t.{column_name}
                    WHERE l.id IS NULL AND t.{column_name} IS NOT NULL
                """

        cursor.execute(sql)
        data = cursor.fetchone()[0]
        if data != 0:
            sys.exit(f"Found {data} unlinked entries in table {table} after locus_id foreign key check\n")

    db.close()

def read_from_gtf(working_dir, ensembl_gtf, exclude_biotypes):
    genes_output_file = working_dir+"/ensembl_genes_grch38.txt"
    error_log_file = working_dir+"/ensembl_genes_grch38_error.log"

    gene_symbol_2_stable_id = {}
    unique_gene_symbol_2_stable_id = {}
    unique_stable_id_2_gene_symbol = {}

    with gzip.open(ensembl_gtf, "rt", encoding='utf-8') as fh:
        for line in fh:
            if line.startswith("#"):
                continue

            data = line.strip().split("\t")
            attribs_list = {}
            for attrib in data[8].split(";"):
                if attrib.strip().startswith("gene_id") or attrib.strip().startswith("gene_name") or attrib.strip().startswith("gene_biotype"):
                    get_attrib_values = attrib.strip().replace("\"", "").split(" ")
                    if get_attrib_values and len(get_attrib_values) == 2:
                        attribs_list[get_attrib_values[0]] = get_attrib_values[1]

            if ("gene_id" in attribs_list and attribs_list["gene_id"] and
                "gene_name" in attribs_list and attribs_list["gene_name"] and
                "gene_biotype" in attribs_list and attribs_list["gene_biotype"] and
                not any(re.search(sub, attribs_list["gene_biotype"]) for sub in exclude_biotypes)):

                if attribs_list["gene_name"] not in gene_symbol_2_stable_id:
                    gene_symbol_2_stable_id[attribs_list["gene_name"]] = set()
                    gene_symbol_2_stable_id[attribs_list["gene_name"]].add(attribs_list["gene_id"])
                else:
                    gene_symbol_2_stable_id[attribs_list["gene_name"]].add(attribs_list["gene_id"])

    # print("-> size:", len(gene_symbol_2_stable_id))

    with open(error_log_file, "w") as wr, open(genes_output_file, "w") as wr_genes:
        for gene_symbol in gene_symbol_2_stable_id:
            if len(gene_symbol_2_stable_id[gene_symbol]) > 1:
                wr.write(f"ERROR: more than one stable_id for {gene_symbol}: {gene_symbol_2_stable_id[gene_symbol]}\n")
            else:
                unique_gene_symbol_2_stable_id[gene_symbol] = list(gene_symbol_2_stable_id[gene_symbol])[0]
                wr_genes.write(f"{unique_gene_symbol_2_stable_id[gene_symbol]}\t{gene_symbol}\n")
                # save the ensembl id as key for easier comparison later
                unique_stable_id_2_gene_symbol[unique_gene_symbol_2_stable_id[gene_symbol]] = gene_symbol

    return unique_gene_symbol_2_stable_id, unique_stable_id_2_gene_symbol

def get_g2p_genes(db_host, db_port, db_name, user, password):
    g2p_genes = {}
    g2p_genes_not_used = {}

    sql =   """
                SELECT l.name, li.identifier FROM locus l
                LEFT JOIN locus_identifier li ON li.locus_id = l.id
                LEFT JOIN source s ON s.id = li.source_id
                WHERE s.name = 'Ensembl'
            """

    sql_not_used =   """
                        SELECT l.name, li.identifier FROM locus l
                        LEFT JOIN locus_identifier li ON l.id = li.locus_id
                        LEFT JOIN source s ON s.id = li.source_id
                        LEFT JOIN locus_genotype_disease lgd ON lgd.locus_id = l.id
                        LEFT JOIN gene_stats g ON g.gene_id = l.id
                        LEFT JOIN uniprot_annotation u ON u.gene_id = l.id
                        LEFT JOIN gene_disease d ON d.gene_id = l.id
                        WHERE l.id is not null and lgd.locus_id is null and g.gene_id is null and 
                        u.gene_id is null and d.gene_id is null and s.name = 'Ensembl'
                    """

    db = MySQLdb.connect(host=db_host, port=db_port, user=user, passwd=password, db=db_name)
    cursor = db.cursor()
    cursor.execute(sql)
    data = cursor.fetchall()
    for row in data:
        g2p_genes[row[1]] = row[0]
    
    cursor.execute(sql_not_used)
    data_other = cursor.fetchall()
    for row in data_other:
        g2p_genes_not_used[row[1]] = row[0]

    db.close()

    return g2p_genes, g2p_genes_not_used


def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--config", required=True, help="Config file")
    parser.add_argument("--version", required=True, help="Ensembl release version")
    parser.add_argument("--working_dir", required=True, help="Working directory")
    args = parser.parse_args()

    config_file = args.config
    version = args.version
    working_dir = args.working_dir

    hgnc_file_url = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
    hgnc_file = working_dir+"/hgnc_complete_set.txt"
    ensembl_gtf_url = f"https://ftp.ensembl.org/pub/release-{version}/gtf/homo_sapiens/Homo_sapiens.GRCh38.{version}.chr.gtf.gz"
    ensembl_gtf = working_dir+f"/Homo_sapiens.GRCh38.{version}.chr.gtf.gz"

    if not os.path.exists(working_dir):
        sys.exit(f"Invalid directory '{working_dir}'")

    try:
        urllib.request.urlretrieve(hgnc_file_url, hgnc_file)
    except:
        sys.exit(f"Problem while fetching HGNC file '{hgnc_file_url}'")

    try:
        urllib.request.urlretrieve(ensembl_gtf_url, ensembl_gtf)
    except:
        sys.exit(f"Problem while fetching Ensembl GTF file '{ensembl_gtf_url}'")

    # Load the config file
    config = configparser.ConfigParser()
    config.read(config_file)

    try:
        g2p_db_details = config['g2p_database']
    except KeyError:
        sys.exit("Config: G2P database details are missing from the config file")
    else:
        db_host = g2p_db_details['host']
        db_port = int(g2p_db_details['port'])
        db_name = g2p_db_details['database']
        user = g2p_db_details['user']
        password = g2p_db_details['password']

    exclude_biotypes = [
        "pseudogene",
        "misc_RNA"
    ]

    g2p_tables_with_locus_id_link = [
        "locus_identifier",
        "locus_attrib",
        "locus_genotype_disease",
        "uniprot_annotation",
        "gene_stats",
        "gene_disease",
        "gene2phenotype_app_historicallocusgenotypedisease" # history table
    ]

    # Check each table with locus id if all data in that column exists in the locus table
    # locus_id_foreign_key_check(db_host, db_port, db_name, user, password, g2p_tables_with_locus_id_link)

    # Get the current G2P genes (gene symbol and corresponding Ensembl ID)
    get_g2p_gene_ids, g2p_genes_not_used = get_g2p_genes(db_host, db_port, db_name, user, password)

    # Get the Ensembl genes from the gtf file
    # Only considers genes that are linked to one Ensembl ID
    unique_gene_symbol_2_stable_id, unique_stable_id_2_gene_symbol = read_from_gtf(working_dir, ensembl_gtf, exclude_biotypes)

    # Check if there are G2P genes that should be removed because they are not valid in the gtf file
    for stable_id in get_g2p_gene_ids:
        if stable_id not in unique_stable_id_2_gene_symbol:
            gene_symbol = get_g2p_gene_ids[stable_id]
            print(f"---{stable_id} ({gene_symbol}) not found in gtf")
            if stable_id in g2p_genes_not_used:
                print(f"-> can deleted {stable_id} ({gene_symbol}) from g2p")

if __name__ == '__main__':
    main()