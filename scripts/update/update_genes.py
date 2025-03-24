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
    gene_symbol_details = {}
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
                    # Save the chr details for each gene symbol
                    gene_symbol_details[attribs_list["gene_name"]] = {
                        "chr": data[0],
                        "start": data[3],
                        "end": data[4],
                        "strand": data[6]
                    }
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
                unique_stable_id_2_gene_symbol[unique_gene_symbol_2_stable_id[gene_symbol]] = {
                    "gene_symbol": gene_symbol,
                    "chr": gene_symbol_details[gene_symbol]["chr"],
                    "start": gene_symbol_details[gene_symbol]["start"],
                    "end": gene_symbol_details[gene_symbol]["end"],
                    "strand": gene_symbol_details[gene_symbol]["strand"]
                }

    return unique_gene_symbol_2_stable_id, unique_stable_id_2_gene_symbol

def get_g2p_genes(db_host, db_port, db_name, user, password):
    g2p_genes = {}
    g2p_genes_by_symbol = {}
    g2p_genes_not_used = {}

    sql =   """
                SELECT l.name, li.identifier, l.id, l.start, l.end, seq.name FROM locus l
                LEFT JOIN locus_identifier li ON li.locus_id = l.id
                LEFT JOIN source s ON s.id = li.source_id
                LEFT JOIN sequence seq ON seq.id = l.sequence_id
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
        g2p_genes[row[1]] = {
            "gene_symbol": row[0],
            "locus_id": row[2],
            "locus_start": row[3],
            "locus_end": row[4],
            "locus_chr": row[5],
        }
        g2p_genes_by_symbol[row[0]] = row[1]

    cursor.execute(sql_not_used)
    data_other = cursor.fetchall()
    for row in data_other:
        g2p_genes_not_used[row[1]] = row[0]

    db.close()

    return g2p_genes, g2p_genes_by_symbol, g2p_genes_not_used

def update_genes(g2p_gene_ids, g2p_genes_by_symbol, unique_stable_id_2_gene_symbol, db_host, db_port, db_name, user, password, working_dir):
    update_report = working_dir+"/report_gene_updates.txt"
    new_genes_report = working_dir+"/report_new_genes.txt"

    sequence_chrs = {}

    sql_sel_sequences = """
                            SELECT id, name FROM sequence
                        """

    sql_sel_attribs =   """
                            SELECT id FROM attrib_type WHERE code = 'gene_synonym'
                        """
    
    sql_sel_attrib_gene =   """
                                SELECT id FROM attrib WHERE value = 'gene'
                            """
    
    sql_sel_source =    """
                            SELECT id FROM source WHERE name = 'Ensembl'
                        """

    sql_update_gene =   """
                            UPDATE locus SET name = %s WHERE id = %s
                        """

    sql_insert_attrib = """
                            INSERT INTO locus_attrib(value, is_deleted, attrib_type_id, locus_id, source_id)
                            VALUES(%s, %s, %s, %s, %s)
                        """

    sql_update_stable_id =  """
                                UPDATE locus_identifier SET identifier = %s
                                WHERE source_id = %s AND locus_id = %s
                            """

    sql_insert_gene =   """
                            INSERT INTO locus(start, end, strand, name, type_id, sequence_id)
                            VALUES(%s, %s, %s, %s, %s, %s)
                        """

    db = MySQLdb.connect(host=db_host, port=db_port, user=user, passwd=password, db=db_name)
    cursor = db.cursor()
    # Get all sequences (chrs)
    cursor.execute(sql_sel_sequences)
    data = cursor.fetchall()
    for row in data:
        sequence_chrs[row[1]] = row[0]

    cursor.execute(sql_sel_attribs)
    attrib_type_id = cursor.fetchone()
    cursor.execute(sql_sel_attrib_gene)
    gene_attrib_id = cursor.fetchone()
    cursor.execute(sql_sel_source)
    source_id = cursor.fetchone()

    with open(update_report, "w") as wr, open(new_genes_report, "w") as wr_new:
        for stable_id in unique_stable_id_2_gene_symbol:
            new_gene_symbol = unique_stable_id_2_gene_symbol[stable_id]["gene_symbol"]

            try:
                g2p_gene_symbol = g2p_gene_ids[stable_id]["gene_symbol"]
            except:
                # The stable_id is not found in G2P, it could mean one of the following:
                # 1) gene is in g2p but in the gtf file the stable id has been updated for the gene symbol
                # 2) this is a new gene
                if new_gene_symbol in g2p_genes_by_symbol:
                    # Scenario 1
                    g2p_current_stable_id = g2p_genes_by_symbol[new_gene_symbol]
                    cursor.execute(sql_update_stable_id, [stable_id, source_id, g2p_gene_ids[g2p_current_stable_id]["locus_id"]])
                    db.commit()
                    # Write to report
                    wr_new.write(f"UPDATE: locus_id = {g2p_gene_ids[g2p_current_stable_id]['locus_id']} gene symbol {new_gene_symbol} new stable_id {stable_id} (previous stable_id {g2p_genes_by_symbol[new_gene_symbol]})\n")
                else:
                    # Scenario 2
                    if unique_stable_id_2_gene_symbol[stable_id]["strand"] == "-":
                        strand = -1
                    else:
                        strand = 1

                    cursor.execute(sql_insert_gene, [
                        unique_stable_id_2_gene_symbol[stable_id]["start"],
                        unique_stable_id_2_gene_symbol[stable_id]["end"],
                        strand,
                        new_gene_symbol,
                        gene_attrib_id,
                        sequence_chrs[unique_stable_id_2_gene_symbol[stable_id]["chr"]]
                    ])
                    db.commit()
                    # Write to report
                    wr_new.write(f"ADD: new gene {new_gene_symbol} stable_id {stable_id}\n")
            else:
                # The stable_id is found in G2P associated with a different gene symbol
                if g2p_gene_symbol != new_gene_symbol:
                    # Check if we can update the gene symbol
                    if ((looks_like_identifier(new_gene_symbol) and not looks_like_identifier(g2p_gene_symbol)) or
                        new_gene_symbol in g2p_gene_ids):
                        # Add gene symbol as a synonym (locus_attrib)
                        cursor.execute(sql_insert_attrib, [new_gene_symbol, 0, attrib_type_id, g2p_gene_ids[stable_id]["locus_id"], source_id])
                        db.commit()
                        # Write to report
                        wr.write(f"ADD SYNONYM: locus_id = {g2p_gene_ids[stable_id]['locus_id']} add locus_attrib {new_gene_symbol}\n\n")
                    else:
                        # Update gene symbol
                        cursor.execute(sql_update_gene, [new_gene_symbol, g2p_gene_ids[stable_id]["locus_id"]])
                        # Add gene symbol as a synonym (locus_attrib)
                        cursor.execute(sql_insert_attrib, [g2p_gene_symbol, 0, attrib_type_id, g2p_gene_ids[stable_id]["locus_id"], source_id])
                        db.commit()
                        # Write to report
                        wr.write(f"UPDATE: locus_id = {g2p_gene_ids[stable_id]['locus_id']} gene symbol {g2p_gene_symbol} updated to {new_gene_symbol}\n")
                        wr.write(f"ADD SYNONYM: locus_id = {g2p_gene_ids[stable_id]['locus_id']} add locus_attrib {g2p_gene_symbol}\n\n")

    db.close()

def looks_like_identifier(symbol):
    return bool(re.match(r'^[A-Z]+[0-9]+\.[0-9]+', symbol))


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
    locus_id_foreign_key_check(db_host, db_port, db_name, user, password, g2p_tables_with_locus_id_link)

    # Get the current G2P genes (Ensembl ID and corresponding gene symbol) -> g2p_gene_ids
    # Get the current G2P genes not used in any table -> g2p_genes_not_used
    g2p_gene_ids, g2p_genes_by_symbol, g2p_genes_not_used = get_g2p_genes(db_host, db_port, db_name, user, password)

    # Get the Ensembl genes from the gtf file
    # Only considers genes that are linked to one Ensembl ID
    # ensembl_gtf = working_dir+f"/Homo_sapiens.GRCh38.114.chr.gtf.gz"
    unique_gene_symbol_2_stable_id, unique_stable_id_2_gene_symbol = read_from_gtf(working_dir, ensembl_gtf, exclude_biotypes)

    # Update gene_symbol or if gene not in G2P insert new gene
    # don't update gene_symbol where readable gene_symbol has been replaced by e.g. AC080038.1 TODO: and gene_symbol is used by G2P
    # If the gene_symbol is already present in the locus table only update the ensembl stable_id
    update_genes(g2p_gene_ids, g2p_genes_by_symbol, unique_stable_id_2_gene_symbol, db_host, db_port, db_name, user, password, working_dir)

    # TODO
    # # Check if there are G2P genes that should be removed because they are not valid in the gtf file
    # for stable_id in g2p_gene_ids:
    #     if stable_id not in unique_stable_id_2_gene_symbol:
    #         gene_symbol = g2p_gene_ids[stable_id]["gene_symbol"]
    #         print(f"---{stable_id} ({gene_symbol}) not found in gtf")
    #         if stable_id in g2p_genes_not_used:
    #             print(f"-> can deleted {stable_id} ({gene_symbol}) from g2p")

if __name__ == '__main__':
    main()