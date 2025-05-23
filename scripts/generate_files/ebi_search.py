#!/usr/bin/env python3

import sys
import argparse
import MySQLdb
import configparser
from datetime import datetime
from lxml import etree as ET


def get_g2p_version(db_host: str, db_port: int, db_name: str, user: str, password: str) -> str:
    """
    Fetch the latest G2P version from the database.

    Args:
        db_host (str): G2P database host name
        db_port (int): G2P database port number
        db_name (str): G2P database name
        user (str): user with read access
        password (str): password

    Returns:
        str: G2P version
    """    
    try:
        db = MySQLdb.connect(host=db_host, port=db_port, user=user, passwd=password, db=db_name)
        cursor = db.cursor()

        sql = """   SELECT m.version
                    FROM meta m
                    LEFT JOIN source s ON s.id = m.source_id
                    WHERE s.name = 'G2P'
                    ORDER BY m.date_update DESC LIMIT 1
            """

        cursor.execute(sql)
        g2p_version = cursor.fetchone()
        cursor.close()
        db.close()

    except MySQLdb.Error as error:
        print("Database connection failed:", error)
        sys.exit(1)

    return g2p_version[0]


def dump_g2p_records(db_host: str, db_port: int, db_name: str, user: str, password: str) -> dict[str, dict]:
    """
    Queries the G2P database to dump all records associated with public panels.

    Args:
        db_host (str): G2P database host name
        db_port (int): G2P database port number
        db_name (str): G2P database name
        user (str): user with read access
        password (str): password

    Returns:
        dict[str, dict]: All records associated with public panels
    """
    records = {}

    try:
        db = MySQLdb.connect(host=db_host, port=db_port, user=user, passwd=password, db=db_name)
        cursor = db.cursor()

        # Select distinct records
        sql = """   SELECT g2p.stable_id, d.name, a1.value, a2.value, m.value, l.name,
                    GROUP_CONCAT(DISTINCT p.name SEPARATOR ',') AS panel_names,
                    GROUP_CONCAT(DISTINCT li.identifier SEPARATOR ',') AS identifiers,
                    GROUP_CONCAT(DISTINCT o.accession SEPARATOR ',') AS ontology_terms
                    FROM locus_genotype_disease lgd
                    LEFT JOIN g2p_stableid g2p ON g2p.id = lgd.stable_id
                    LEFT JOIN lgd_panel panel ON panel.lgd_id = lgd.id
                    LEFT JOIN panel p ON p.id = panel.panel_id
                    LEFT JOIN disease d ON d.id = lgd.disease_id
                    LEFT JOIN disease_ontology_term do ON do.disease_id = d.id
                    LEFT JOIN ontology_term o ON o.id = do.ontology_term_id
                    LEFT JOIN attrib a1 ON a1.id = lgd.confidence_id
                    LEFT JOIN attrib a2 ON a2.id = lgd.genotype_id
                    LEFT JOIN cv_molecular_mechanism m ON m.id = lgd.mechanism_id
                    LEFT JOIN locus l ON l.id = lgd.locus_id
                    LEFT JOIN locus_identifier li ON li.locus_id = l.id
                    WHERE p.is_visible = 1
                    GROUP BY g2p.stable_id, d.name, a1.value, a2.value, m.value, l.name
            """

        # TODO: add publications
        # cross reference: EUROPEPMC
        # https://europepmc.org/article/MED/{pmid}

        cursor.execute(sql)
        data = cursor.fetchall()
        for row in data:
            g2p_id = row[0]
            panels = row[6].split(",")
            gene_ids = row[7].split(",")

            if g2p_id not in records:
                records[g2p_id] = {
                    "disease": row[1],
                    "confidence": row[2],
                    "genotype": row[3],
                    "mechanism": row[4],
                    "gene": row[5],
                    "panels": panels,
                    "gene_ids": gene_ids
                }
                # Check if there are ontology terms
                if row[8]:
                    records[g2p_id]["ontology_terms"] = row[8].split(",")

            else:
                print(f"WARNING: duplicated record '{g2p_id}'")

        cursor.close()
        db.close()

    except MySQLdb.Error as error:
        print("Database connection failed:", error)
        sys.exit(1)

    return records


def add_elem(parent: object, name: str, value: str) -> None:
    """
    Add a field to the XML obj.

    Args:
        parent (object): ET element obj
        name (str): name of the element field
        value (str): value of the field
    """
    if value:
        child = ET.SubElement(parent, name)
        child.text = value


def create_xml(g2p_version: str, g2p_records: dict[str, dict]) -> bytes:
    """
    Generate the XML data with all the records that we want to display on the
    EBI search engine.

    Args:
        g2p_version (str): G2P version
        g2p_records (dict[str, dict]): List of all records linked to public panels

    Returns:
        bytes: xml data generated with the G2P records info
    """    
    database_elem = ET.Element("database")

    # Add metadata
    add_elem(database_elem, "name", "G2P")
    add_elem(database_elem, "description", "Gene2Phenotype (G2P) is a detailed collection of expert curated gene-disease associations with information on allelic requirement, observed variant classes and disease mechanism")
    add_elem(database_elem, "url", "https://www.ebi.ac.uk/gene2phenotype/")
    add_elem(database_elem, "url_search", "https://www.ebi.ac.uk/gene2phenotype/lgd/")
    add_elem(database_elem, "release", g2p_version)
    add_elem(database_elem, "entry_count", str(len(g2p_records)))

    # Add records
    entries_elem = ET.SubElement(database_elem, "entries")
    for entry in g2p_records:
        entry_elem = ET.SubElement(entries_elem, "entry", id=entry, acc=entry)
        entry_name = f"{g2p_records[entry]['disease']} ({g2p_records[entry]['confidence']})"
        add_elem(entry_elem, "name", entry_name)
        # Additional fields
        add_fields_elem = ET.SubElement(entry_elem, "additional_fields")
        f = ET.SubElement(add_fields_elem, "field", name="gene")
        f.text = g2p_records[entry]["gene"]
        f = ET.SubElement(add_fields_elem, "field", name="disease")
        f.text = g2p_records[entry]["disease"]
        f = ET.SubElement(add_fields_elem, "field", name="genotype")
        f.text = g2p_records[entry]["genotype"]
        f = ET.SubElement(add_fields_elem, "field", name="mechanism")
        f.text = g2p_records[entry]["mechanism"]
        f = ET.SubElement(add_fields_elem, "field", name="confidence")
        f.text = g2p_records[entry]["confidence"]
        # Cross references - gene ID
        xrefs_elem = ET.SubElement(entry_elem, "cross_references")
        for xref in g2p_records[entry]["gene_ids"]:
            if xref.startswith("ENSG"):
                db = "ENSEMBL_GENE"
            elif xref.startswith("HGNC:"):
                db = "HGNC"
            else:
                # Skip the OMIM gene ID because EBI search does not support this ID
                continue
            ET.SubElement(xrefs_elem, "ref", dbname=db, dbkey=xref)
        # Cross references - disease ontology
        try:
            ontology_list = g2p_records[entry]["ontology_terms"]
        except:
            continue
        else:
            for ontology_xref in ontology_list:
                if ontology_xref.startswith("MONDO"):
                    db = "Mondo"
                elif ontology_xref.startswith("Orphanet"):
                    db = "Orphanet"
                else:
                    db = "OMIM_DISEASE"
                ET.SubElement(xrefs_elem, "ref", dbname=db, dbkey=ontology_xref)

    return ET.tostring(database_elem, pretty_print=True, encoding="UTF-8", xml_declaration=True)


def main():
    parser = argparse.ArgumentParser(description="Generate a G2P XML file for the EBI search engine")
    parser.add_argument("--config", required=True, help="Config file containing details to the G2P database")
    args = parser.parse_args()

    config_file = args.config

    # Load the config file
    config = configparser.ConfigParser()
    config.read(config_file)

    try:
        g2p_config = config['g2p_database']
    except KeyError:
        sys.exit("ERROR: 'g2p_database' missing from config file")
    else:
        db_host = g2p_config['host']
        db_port = g2p_config['port']
        db_name = g2p_config['name']
        user = g2p_config['user']
        password = g2p_config['password']

    # Get the G2P version from the meta table
    g2p_version = get_g2p_version(db_host, int(db_port), db_name, user, password)
    # Fetch all records to be available in the EBI search
    g2p_records = dump_g2p_records(db_host, int(db_port), db_name, user, password)
    # Generate the XML with all records
    xml_data = create_xml(g2p_version, g2p_records)

    # Write to the output file - open in binary mode because 'xml_data' is a bytes obj
    output_file = f"g2p_records.xml"

    with open(output_file, "wb") as wr:
        wr.write(xml_data)


if __name__ == '__main__':
    main()