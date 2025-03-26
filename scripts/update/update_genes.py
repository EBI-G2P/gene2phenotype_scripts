#!/usr/bin/env python3

import os.path
import sys
import argparse
import MySQLdb
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

def read_from_gtf(ensembl_gtf, exclude_biotypes):
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

    return unique_stable_id_2_gene_symbol

def get_g2p_genes(db_host, db_port, db_name, user, password):
    g2p_genes = {}
    g2p_genes_by_symbol = {}

    sql =   """
                SELECT l.name, li.identifier, l.id, l.start, l.end, seq.name FROM locus l
                LEFT JOIN locus_identifier li ON li.locus_id = l.id
                LEFT JOIN source s ON s.id = li.source_id
                LEFT JOIN sequence seq ON seq.id = l.sequence_id
                WHERE s.name = 'Ensembl'
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

    db.close()

    return g2p_genes, g2p_genes_by_symbol

def get_g2p_genes_hgnc(db_host, db_port, db_name, user, password):
    """
        Dump the full G2P gene set.
        Returns a dictionary where the key is the gene symbol.
        Called by: update_xrefs()
    """
    g2p_genes_by_symbol = {}

    sql =   """
                SELECT l.name, li.identifier, l.id, la.value FROM locus l
                LEFT JOIN locus_identifier li ON li.locus_id = l.id
                LEFT JOIN locus_attrib la ON la.locus_id = l.id
                LEFT JOIN source s ON s.id = li.source_id
                WHERE s.name = 'Ensembl'
            """

    sql_ids =   """
                    SELECT l.name, li.identifier FROM locus l
                    LEFT JOIN locus_identifier li ON li.locus_id = l.id
                """

    db = MySQLdb.connect(host=db_host, port=db_port, user=user, passwd=password, db=db_name)
    cursor = db.cursor()
    cursor.execute(sql)
    data = cursor.fetchall()
    for row in data:
        if row[0] not in g2p_genes_by_symbol:
            g2p_genes_by_symbol[row[0]] = {
                "stable_id": row[1],
                "locus_id": row[2],
                "synonyms": [row[3]]
            }
        else:
            g2p_genes_by_symbol[row[0]]["synonyms"].append(row[3])

    cursor.execute(sql_ids)
    data_ids = cursor.fetchall()
    for row in data_ids:
        if row[1].startswith("HGNC"):
            g2p_genes_by_symbol[row[0]]["hgnc_id"] = row[1]
        elif row[1].isdigit():
            g2p_genes_by_symbol[row[0]]["omim_id"] = row[1]

    db.close()

    return g2p_genes_by_symbol

def get_g2p_genes_not_used(db_host, db_port, db_name, user, password):
    g2p_genes_not_used = {}

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
    cursor.execute(sql_not_used)
    data_other = cursor.fetchall()
    for row in data_other:
        g2p_genes_not_used[row[1]] = row[0]

    db.close()

    return g2p_genes_not_used

def update_genes(g2p_gene_ids, g2p_genes_by_symbol, unique_stable_id_2_gene_symbol, db_host, db_port, db_name, user, password):
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

    sql_update_coord =  """
                            UPDATE locus SET start = %s, end = %s, sequence_id = %s WHERE id = %s
                        """

    sql_insert_attrib = """
                            INSERT INTO locus_attrib(value, is_deleted, attrib_type_id, locus_id, source_id)
                            VALUES(%s, %s, %s, %s, %s)
                        """

    sql_insert_identifier = """
                                INSERT INTO locus_identifier(identifier, description, locus_id, source_id)
                                VALUES(%s, %s, %s, %s)
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

                    # Add gene
                    cursor.execute(sql_insert_gene, [
                        unique_stable_id_2_gene_symbol[stable_id]["start"],
                        unique_stable_id_2_gene_symbol[stable_id]["end"],
                        strand,
                        new_gene_symbol,
                        gene_attrib_id,
                        sequence_chrs[unique_stable_id_2_gene_symbol[stable_id]["chr"]]
                    ])
                    db.commit()
                    locus_id = cursor.lastrowid
                    # Add Ensembl stable_id
                    cursor.execute(sql_insert_identifier, [stable_id, None, locus_id, source_id])
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
                        # Add old gene symbol as a synonym (locus_attrib)
                        cursor.execute(sql_insert_attrib, [g2p_gene_symbol, 0, attrib_type_id, g2p_gene_ids[stable_id]["locus_id"], source_id])
                        db.commit()
                        # Write to report
                        wr.write(f"\nUPDATE: locus_id = {g2p_gene_ids[stable_id]['locus_id']} gene symbol {g2p_gene_symbol} updated to {new_gene_symbol}\n")
                        wr.write(f"ADD SYNONYM: locus_id = {g2p_gene_ids[stable_id]['locus_id']} add locus_attrib {g2p_gene_symbol}\n")
                    
                    # Compare the coordinates
                    if (g2p_gene_ids[stable_id]["start"] != unique_stable_id_2_gene_symbol[stable_id]["start"] or
                        g2p_gene_ids[stable_id]["end"] != unique_stable_id_2_gene_symbol[stable_id]["end"] or
                        g2p_gene_ids[stable_id]["chr"] != unique_stable_id_2_gene_symbol[stable_id]["chr"]):
                        cursor.execute(sql_update_coord, [
                            unique_stable_id_2_gene_symbol[stable_id]["start"],
                            unique_stable_id_2_gene_symbol[stable_id]["end"],
                            sequence_chrs[unique_stable_id_2_gene_symbol[stable_id]["chr"]],
                            g2p_gene_ids[stable_id]['locus_id']
                        ])
                        db.commit()
                        wr.write(f"UPDATE COORD: locus_id = {g2p_gene_ids[stable_id]['locus_id']} (previous {g2p_gene_ids[stable_id]['chr']}:{g2p_gene_ids[stable_id]['start']}-{g2p_gene_ids[stable_id]['end']})\n")

    db.close()

def looks_like_identifier(symbol):
    return bool(re.match(r'^[A-Z]+[0-9]+\.[0-9]+', symbol))

def update_xrefs(hgnc_file, db_host, db_port, db_name, user, password):
    report_hgnc_file = working_dir+"/report_hgnc_updates.txt"

    # Get a new dump of the g2p genes - after they were updated with the gft file
    g2p_genes_by_symbol = get_g2p_genes_hgnc(db_host, db_port, db_name, user, password)

    mappings = {}

    with open(hgnc_file, "r",) as fh:
        for line in fh:
            if line.startswith("hgnc_id"):
                header = line.strip().split("\t")
                continue

            data = line.strip().split("\t")
            data = dict(zip(header, data))
            gene_symbol = data["symbol"]

            for id_type in ["hgnc_id", "prev_symbol", "ensembl_gene_id", "omim_id"]:
                try:
                    value = data[id_type]
                except:
                    continue
                else:
                    # clean the value before saving it
                    value = value.replace("\"", "")
                    if gene_symbol not in mappings:
                        mappings[gene_symbol] = {}
                        mappings[gene_symbol][id_type] = [value]
                    elif id_type in mappings[gene_symbol]:
                        mappings[gene_symbol][id_type].append(value)
                    else:
                        mappings[gene_symbol][id_type] = [value]

    # Update
    with open(report_hgnc_file, "w",) as wr:
        for symbol in mappings:
            # Get current g2p data for this gene symbol
            try: 
                g2p_data = g2p_genes_by_symbol[symbol]
            except:
                continue
            else:
                # Update or add HGNC ID
                if len(mappings[symbol]["hgnc_id"]) == 1:
                    if "hgnc_id" in g2p_data and mappings[symbol]["hgnc_id"][0] != g2p_data["hgnc_id"]:
                        update_g2p_identifier("HGNC_update", mappings[symbol]["hgnc_id"][0], g2p_data["locus_id"], db_host, db_port, db_name, user, password)
                        wr.write(f"UPDATE HGNC: locus_id = {g2p_data['locus_id']} gene symbol {symbol} {mappings[symbol]['hgnc_id'][0]}\n")
                    elif "hgnc_id" not in g2p_data:
                        update_g2p_identifier("HGNC_insert", mappings[symbol]["hgnc_id"][0], g2p_data["locus_id"], db_host, db_port, db_name, user, password)
                        wr.write(f"ADD HGNC: locus_id = {g2p_data['locus_id']} gene symbol {symbol} {mappings[symbol]['hgnc_id'][0]}\n")

                # Update or add OMIM ID
                if "omim_id" in mappings[symbol]:
                    # Format OMIM IDs
                    omim_ids = mappings[symbol]["omim_id"][0].split("|")
                    if len(omim_ids) == 1 and omim_ids[0] != "" and "omim_id" in g2p_data and omim_ids[0] != g2p_data["omim_id"]:
                        update_g2p_identifier("OMIM_update", omim_ids[0], g2p_data["locus_id"], db_host, db_port, db_name, user, password)
                        wr.write(f"UPDATE OMIM ID: locus_id = {g2p_data['locus_id']} gene symbol {symbol} {omim_ids[0]}\n")
                    elif "omim_id" not in g2p_data and omim_ids[0] != "" and len(omim_ids) == 1:
                        update_g2p_identifier("OMIM_insert", omim_ids[0], g2p_data["locus_id"], db_host, db_port, db_name, user, password)
                        wr.write(f"ADD OMIM ID: locus_id = {g2p_data['locus_id']} gene symbol {symbol} {omim_ids[0]}\n")

                # Add gene synonyms (locus_attrib)
                if "prev_symbol" in mappings[symbol]:
                    prev_symbols = mappings[symbol]["prev_symbol"][0].split("|")
                    for gene_syn in prev_symbols:
                        if (gene_syn.upper() not in g2p_data["synonyms"] and gene_syn not in g2p_data["synonyms"]) and gene_syn != "":
                            add_gene_synonym(gene_syn, g2p_data['locus_id'], db_host, db_port, db_name, user, password)
                            wr.write(f"ADD GENE PREV SYMBOL: locus_id = {g2p_data['locus_id']} gene symbol {symbol} prev symbol {gene_syn}\n")

def update_g2p_identifier(type, new_hgnc_id, locus_id, db_host, db_port, db_name, user, password):
    """
        type: HGNC or OMIM
    """
    source_name, query_type = type.split("_")

    sql_source =    """
                        SELECT id FROM source
                        WHERE name = %s
                    """

    sql_update_hgnc =   """
                            UPDATE locus_identifier SET identifier = %s
                            WHERE locus_id = %s AND source_id = %s
                        """

    sql_insert_hgnc =   """
                            INSERT INTO locus_identifier(identifier, description, locus_id, source_id)
                            VALUES(%s, %s, %s, %s)
                        """

    db = MySQLdb.connect(host=db_host, port=db_port, user=user, passwd=password, db=db_name)
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

def add_gene_synonym(symbol, locus_id, db_host, db_port, db_name, user, password):

    sql_sel_attribs =   """
                            SELECT id FROM attrib_type WHERE code = 'gene_synonym'
                        """

    sql_sel_source =    """
                            SELECT id FROM source WHERE name = 'Ensembl'
                        """

    sql_insert_syn =   """
                            INSERT INTO locus_attrib(value, is_deleted, attrib_type_id, locus_id, source_id)
                            VALUES(%s, %s, %s, %s, %s)
                        """

    db = MySQLdb.connect(host=db_host, port=db_port, user=user, passwd=password, db=db_name)
    cursor = db.cursor()
    cursor.execute(sql_sel_attribs)
    attrib_type_id = cursor.fetchone()
    cursor.execute(sql_sel_source)
    source_id = cursor.fetchone()

    # Insert gene synonym
    cursor.execute(sql_insert_syn, [symbol, 0, attrib_type_id, locus_id, source_id])
    db.commit()
    db.close()


def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--config", required=True, help="Config file")
    parser.add_argument("--version", required=True, help="Ensembl release version")
    parser.add_argument("--working_dir", required=True, help="Working directory")
    args = parser.parse_args()

    config_file = args.config
    version = args.version
    global working_dir
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
    # Get the current G2P genes (gene symbol and corresponding Ensembl ID) -> g2p_genes_by_symbol
    g2p_gene_ids, g2p_genes_by_symbol = get_g2p_genes(db_host, db_port, db_name, user, password)

    # Get the Ensembl genes from the gtf file
    # Only considers genes that are linked to one Ensembl ID
    unique_stable_id_2_gene_symbol = read_from_gtf(ensembl_gtf, exclude_biotypes)

    # Update gene_symbol or if gene not in G2P insert new gene
    # don't update gene_symbol where readable gene_symbol has been replaced by e.g. AC080038.1 TODO: and gene_symbol is used by G2P
    # If the gene_symbol is already present in the locus table only update the ensembl stable_id
    update_genes(g2p_gene_ids, g2p_genes_by_symbol, unique_stable_id_2_gene_symbol, db_host, db_port, db_name, user, password)

    # Update the HGNC IDs
    # This update uses the gene symbol to compare the genes
    # We need a new dump of the genes (after all the previous updates)
    update_xrefs(hgnc_file, db_host, db_port, db_name, user, password)

    # TODO
    # # Check if there are G2P genes that should be removed because they are not valid in the gtf file
    # g2p_genes_not_used = get_g2p_genes_not_used()
    # for stable_id in g2p_gene_ids:
    #     if stable_id not in unique_stable_id_2_gene_symbol:
    #         gene_symbol = g2p_gene_ids[stable_id]["gene_symbol"]
    #         print(f"---{stable_id} ({gene_symbol}) not found in gtf")
    #         if stable_id in g2p_genes_not_used:
    #             print(f"-> can deleted {stable_id} ({gene_symbol}) from g2p")

    # Check again each table with locus id if all data in that column exists in the locus table
    locus_id_foreign_key_check(db_host, db_port, db_name, user, password, g2p_tables_with_locus_id_link)

    # TODO: update meta

if __name__ == '__main__':
    main()