#!/usr/bin/env python3

import sys
import re
import argparse
import datetime
import MySQLdb
import xml.etree.ElementTree as ET
import warnings
import configparser
from pathlib import Path
from os import listdir

"""
    Script      : import_gene_disease.py

    Description : Imports/updates gene-disease associations from OMIM (Ensembl) and Mondo
                  The script can be run in two modes:
                    - import : Runs a full import of the data. This mode should only be run on a newly created db.
                               The script runs a check before this mode is run to make sure the db doesn't include 
                               gene-disease associations beforehand.
                    - update (default) : Updates the data if necessary. This mode should be used in the release cycle.

    Options
                --config: Config file with connection details to Ensembl and G2P databases (mandatory)
                    File format is the following:
                            [ensembl_database]
                            host = <>
                            port = <>
                            user = <>
                            password = <>
                            database = <>

                            [g2p_database]
                            g2p_host = <>
                            g2p_port = <>
                            g2p_database = <>
                            g2p_user = <>
                            g2p_password = <>

                --skip_omim: Option to skip OMIM gene-disease data import (default=off)
                --skip_mondo: Option to skip MONDO gene-disease data import (default=off)
                --mondo_file: MONDO file. Supported format: owl. Only necessary when importing mondo (--mondo).
                --full_import: Option to run a full import of the data. By default, the script runs a data update (default=off)
"""

debug = 0 # to help parse the owl file


def get_g2p_meta(db_host, db_port, db_name, user, password):
    """
        Retrieve meta info from G2P db
        Retrieve the latest OMIM and Mondo updates info
    """

    meta_info = {}

    db = MySQLdb.connect(host=db_host, port=db_port, user=user, passwd=password, db=db_name)
    cursor = db.cursor()

    sql = """ SELECT s.name, m.version, m.date_update
              FROM meta m
              LEFT JOIN source s ON s.id = m.source_id
              WHERE m.`key` = 'import_gene_disease_omim' OR m.`key` = 'import_gene_disease_mondo'
              ORDER BY m.date_update desc
          """

    cursor.execute(sql)
    data = cursor.fetchall()
    for row in data:
        source_name = row[0]
        version = row[1]
        date_update = row[2]

        if source_name not in meta_info or date_update > meta_info[source_name]["date_update"]:
            meta_info[source_name] = { 
                "data_version": version,
                "date_update": date_update
            }

    return meta_info


def get_mim_gene_diseases(db_host, db_port, db_name, user, password):
    """
        Retrieve OMIM gene disease data from Ensembl core db
    """

    gene_diseases = {}

    db = MySQLdb.connect(host=db_host, port=db_port, user=user, passwd=password, db=db_name)
    cursor = db.cursor()
    sql = """ SELECT a.value, g.stable_id, x.dbprimary_acc, x.description 
              FROM external_db e, xref x, object_xref o, gene g, gene_attrib a 
              WHERE e.external_db_id = x.external_db_id AND x.xref_id = o.xref_id 
              AND o.ensembl_id = g.gene_id AND e.db_name = 'MIM_MORBID' AND a.gene_id = g.gene_id 
              AND a.attrib_type_id = 4
          """

    cursor.execute(sql)
    data = cursor.fetchall()
    for row in data:
        disease = row[3].split(';')[0]
        if row[1] not in gene_diseases.keys():
            gene_diseases[row[1]] = [{ 'mim_id': row[2],
                                       'disease': disease }]
        else:
            gene_diseases[row[1]].append({ 'mim_id':row[2],
                                           'disease':disease })

    db.close()
    return gene_diseases


def insert_mim_gene_diseases(db_host, db_port, db_name, user, password, gene_diseases, version):
    """Insert OMIM gene disease data into G2P database"""

    sql_delete = """ DELETE FROM gene_disease
                     WHERE source_id = %s
                 """

    sql_gene = """ SELECT l.id FROM locus l
                   LEFT JOIN locus_identifier i on i.locus_id = l.id
                   LEFT JOIN source s on s.id = i.source_id
                   WHERE i.identifier = %s AND s.name = 'Ensembl'
               """

    sql_source = "SELECT id, name FROM source WHERE name = 'OMIM' OR name = 'Mondo' or name = 'Ensembl'"

    sql_insert = """ INSERT INTO gene_disease(gene_id, disease, identifier, source_id)
                      VALUES(%s, %s, %s, %s)
                 """

    sql_meta = """ INSERT INTO meta(`key`, date_update, is_public, description, source_id, version)
                   VALUES(%s,%s,%s,%s,%s,%s)
               """

    db = MySQLdb.connect(host=db_host, port=db_port, user=user, passwd=password, db=db_name)
    cursor = db.cursor()
    # Fetch source id
    source_ids = {}
    cursor.execute(sql_source)
    data = cursor.fetchall()
    for row in data:
        source_ids[row[1]] = row[0]

    # Delete the existing OMIM data first
    cursor.execute(sql_delete, [source_ids["OMIM"]])

    for stable_id, gd_info in gene_diseases.items():
        cursor.execute(sql_gene, [stable_id])
        data = cursor.fetchone()
        gene_id = None
        if data:
            gene_id = data[0]
        if gene_id is not None:
            for info in gd_info:
                cursor.execute(sql_insert, [gene_id, info['disease'], info['mim_id'], source_ids['OMIM']])

    # Insert import info into meta
    cursor.execute(sql_meta, ['import_gene_disease_omim',
                              datetime.datetime.now(),
                              0,
                              'Import OMIM gene disease associations from Ensembl core db',
                              source_ids['Ensembl'],
                              version])

    db.commit()
    db.close()


def get_mondo_gene_diseases(file):
    """
        Retrieve the gene-disease associations from the Mondo file (owl format).
        This method also retrieves all Mondo IDs excluding 'obsolete' - this data will
        be used to populate table 'disease_external'.

        Returns a tuple of dictionaries with the following structure
        - results (dict)
            key: (string) Mondo ID
            value: (dict) { "disease_name": name of disease, "hgnc_id": HGNC ID }
        - all_diseases (dict)
            key: (string) Mondo ID
            value: (string) name of disease
    """

    results = {}
    all_diseases = {}

    # Open the tmp dir where all the input files are saved
    for tmp_file in listdir("tmp/"):
        with open("tmp/"+tmp_file, "r") as fh:
            mondo_ids = []
            mondo_diseases = {}
            start_class = 0

            for event, elem in ET.iterparse(fh, events=("start", "end")):
                is_element = None
                disease_name = None
                disease_synonym = None
                hgnc_id = None

                if event == "start" and elem.tag == "{http://www.w3.org/2002/07/owl#}Class":
                    # Some entries have a class inside a class
                    # We want to skip those - only the first class is relevant
                    if start_class == 1 and debug == 1:
                        print("Start class when another class still open")

                    if start_class == 0:
                        start_class = 1

                        if debug == 1:
                            ET.dump(elem)

                        for i in elem.iter():
                            # Get mondo ID
                            if i.tag == "{http://www.geneontology.org/formats/oboInOwl#}id" and i.text is not None:
                                is_element = i.text

                            # Get disease name
                            if i.tag == "{http://www.w3.org/2000/01/rdf-schema#}label" and i.text is not None:
                                disease_name = i.text

                            # Get disease synonym
                            # This is necessary if for some reason it can't get the name from the label
                            if i.tag == "{http://www.geneontology.org/formats/oboInOwl#}hasExactSynonym" and i.text is not None:
                                disease_synonym = i.text

                            # Get HGNC ID from the class
                            # This is probably not necessary
                            if(i.tag == "{http://www.w3.org/2002/07/owl#}someValuesFrom"
                            and i.attrib is not None and "{http://www.w3.org/1999/02/22-rdf-syntax-ns#}resource" in i.attrib.keys()
                            and "hgnc" in i.attrib["{http://www.w3.org/1999/02/22-rdf-syntax-ns#}resource"]):
                                hgnc_id = i.attrib["{http://www.w3.org/1999/02/22-rdf-syntax-ns#}resource"].split("/")[-1]

                            if is_element is not None and disease_name is not None and not disease_name.startswith("obsolete"):
                                mondo_ids.append(is_element)
                            
                            # This is a special case
                            elif is_element is not None and disease_name is None and disease_synonym is not None:
                                mondo_ids.append(is_element)
                                disease_name = disease_synonym

                # Class is over
                if event == "end" and elem.tag == "{http://www.w3.org/2002/07/owl#}Class":
                    start_class = 0

                # Save the disease name - this is necessary because the HGNC is outside of the Class
                if is_element is not None and disease_name is not None and is_element not in mondo_diseases:
                    mondo_diseases[is_element] = disease_name

                # Get HGNC ID from outside the elem "Class"
                # At this point the is_element is None - we don't have the Mondo ID anymore
                # but we know that the Mondo ID is the previous ID analysed
                if event == "end" and elem.tag == "{http://www.w3.org/2002/07/owl#}Restriction" and is_element is None and len(mondo_ids) >= 1:
                    is_element = mondo_ids[-1]

                    for i in elem.iter():
                        if(i.tag == "{http://www.w3.org/2002/07/owl#}someValuesFrom" 
                            and "{http://www.w3.org/1999/02/22-rdf-syntax-ns#}resource" in i.attrib.keys()
                            and "hgnc" in i.attrib["{http://www.w3.org/1999/02/22-rdf-syntax-ns#}resource"]):
                            hgnc_id = i.attrib["{http://www.w3.org/1999/02/22-rdf-syntax-ns#}resource"].split("/")[-1]
                            disease_name = mondo_diseases[is_element]

                if is_element is not None and disease_name is not None and hgnc_id is not None and is_element not in results:
                    results[is_element] = { "disease": disease_name, "hgnc_id": hgnc_id }
                if is_element is not None and disease_name is not None and is_element not in results and not disease_name.startswith("obsolete"):
                    all_diseases[is_element] = disease_name

            fh.close()

    return results, all_diseases


def insert_mondo_gene_diseases(db_host, db_port, db_name, user, password, gene_diseases, mondo_version):
    """
    Insert Mondo gene disease data into G2P database
    This method runs a full import
    """
    sql_delete = """ DELETE FROM gene_disease
                     WHERE source_id = %s
                 """

    sql_gene = """ SELECT l.id FROM locus l
                   LEFT JOIN locus_identifier i on i.locus_id = l.id
                   LEFT JOIN source s on s.id = i.source_id
                   WHERE i.identifier = %s AND s.name = 'HGNC'
               """

    sql_source = """ SELECT id, name FROM source WHERE name = 'Mondo'
                 """

    sql_insert = """ INSERT INTO gene_disease(gene_id, disease, identifier, source_id)
                     VALUES(%s, %s, %s, %s)
                 """

    sql_meta = """ INSERT INTO meta(`key`, date_update, is_public, description, source_id, version)
                   VALUES(%s,%s,%s,%s,%s,%s)
               """

    db = MySQLdb.connect(host=db_host, port=db_port, user=user, passwd=password, db=db_name)
    cursor = db.cursor()
    # Fetch source id
    source_id = None
    cursor.execute(sql_source)
    data = cursor.fetchall()
    for row in data:
        source_id = row[0]

    # Delete Mondo data
    cursor.execute(sql_delete, [source_id])

    for gd, gd_info in gene_diseases.items():
        cursor.execute(sql_gene, [f"HGNC:{gd}"])
        data = cursor.fetchall()
        gene_id = None
        for row in data:
            gene_id = row[0]
        if gene_id is not None:
            for info in gd_info:
                cursor.execute(sql_insert, [gene_id, info['disease'], info['mondo_id'], source_id])

    # Insert import info into meta
    cursor.execute(sql_meta, ['import_gene_disease_mondo',
                              datetime.datetime.now(),
                              0,
                              'Import Mondo gene disease associations',
                              source_id,
                              mondo_version])

    db.commit()
    db.close()


def populate_mondo_all_diseases(db_host, db_port, db_name, user, password, mondo_all_diseases):
    sql_ins_mondo = """ INSERT INTO disease_external(disease, identifier, source_id)
                        VALUES(%s, %s, %s)
                    """

    sql_truncate =  """ TRUNCATE TABLE disease_external """

    sql_source = """ SELECT id FROM source WHERE name = 'Mondo' """

    db = MySQLdb.connect(host=db_host, port=db_port, user=user, passwd=password, db=db_name)
    cursor = db.cursor()
    # Fetch source id
    cursor.execute(sql_source)
    data_source = cursor.fetchone()
    source_id = data_source[0]

    # Truncate table
    cursor.execute(sql_truncate)
    db.commit()

    # Prepare input data
    input_data = []
    for mondo_disease in mondo_all_diseases:
        new_data = (mondo_all_diseases[mondo_disease], mondo_disease, source_id)
        input_data.append(new_data)

    # Insert data
    cursor.executemany(sql_ins_mondo, input_data)
    db.commit()
    db.close()


def fetch_mondo_version(file):
    """
        Fetch the Mondo data version from the input file
    """

    version = None

    with open(file, "r") as fh:
            for event, elem in ET.iterparse(fh, events=("start", "end")):
                if event == "start" and elem.tag == "{http://www.w3.org/1999/02/22-rdf-syntax-ns#}RDF":
                    for i in elem.iter():
                        if i.tag == "{http://www.w3.org/2002/07/owl#}versionIRI":
                            resource = i.get("{http://www.w3.org/1999/02/22-rdf-syntax-ns#}resource")
                            version = re.search("[0-9]+\\-[0-9]+\\-[0-9]+", resource)
                            return version.group()

    return version


def split_mondo_file(file):
    Path("tmp").mkdir(parents=True, exist_ok=True)
    output_file = "tmp/mondo_"
    count = 1
    header = None

    with open(file, 'r') as f2:
        data = f2.read()

        for chunk in re.split(r"\s*\<\!\-\-*", data):
            if count == 1:
                header = chunk
                count = count + 1
            else:
                with open(output_file+str(count)+".owl", "w") as wr:
                    wr.write(header)
                    wr.write("\n<!-- "+chunk)
                    if "</rdf:RDF>" not in chunk:
                        wr.write("\n</rdf:RDF>")
                    count = count + 1


def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--config", required=True, help="Config file with details to Ensembl and G2P databases")
    parser.add_argument("--skip_omim", dest="omim", action="store_false", help="Skip importing OMIM gene-disease data (default: off)")
    parser.add_argument("--skip_mondo", dest="mondo", action="store_false", help="Skip importing Mondo gene-disease data (default: off)")
    parser.add_argument("--mondo_file", default='', help="MONDO file in the format: owl")

    args = parser.parse_args()

    config_file = args.config
    mondo_file = args.mondo_file

    # Load the config file
    config = configparser.ConfigParser()
    config.read(config_file)

    db_host = config['ensembl_database']['host']
    db_port = config['ensembl_database']['port']
    db_name = config['ensembl_database']['database']
    user = config['ensembl_database']['user']
    password = config['ensembl_database']['password']
    g2p_db_host = config['g2p_database']['g2p_host']
    g2p_db_port = int(config['g2p_database']['g2p_port'])
    g2p_db_name = config['g2p_database']['g2p_database']
    g2p_user = config['g2p_database']['g2p_user']
    g2p_password = config['g2p_database']['g2p_password']

    # Fetch gene-disease import info from G2P
    g2p_meta_info = get_g2p_meta(g2p_db_host, g2p_db_port, g2p_db_name, g2p_user, g2p_password)

    ### Import or update OMIM ###
    if args.omim:
        # Get version from Ensembl db name
        version = re.search("[0-9]+", db_name)
        print(f"INFO: going to import OMIM data from Ensembl {version.group()}")

        # Compare G2P OMIM last update with Ensembl (OMIM) data version
        if "Ensembl" in g2p_meta_info:
            print(f"INFO: current OMIM data is from Ensembl {g2p_meta_info['Ensembl']['data_version']}")

            if g2p_meta_info["Ensembl"]["data_version"] >= version.group():
                sys.exit("Skipping update: OMIM gene-disease data is up-to-date.")

        # Only fetch OMIM data from Ensembl if import/update can be done
        # Retrive OMIM gene-disease from Ensembl core db
        print("Getting OMIM gene-disease associations...")
        gene_diseases = get_mim_gene_diseases(db_host, int(db_port), db_name, user, password)
        print("Getting OMIM gene-disease associations... done")

        # Insert OMIM gene-disease data into G2P db
        print("Inserting OMIM gene-disease associations into G2P...")
        insert_mim_gene_diseases(g2p_db_host, g2p_db_port, g2p_db_name, g2p_user, g2p_password, gene_diseases, version.group())
        print("Inserting OMIM gene-disease associations into G2P... done")

    ### Import or update Mondo ###
    if args.mondo:
        """
            Input file: it is necessary to perform some cleaning before the import.
            The Mondo file should only have the following type of entries:
                <!-- http://purl.obolibrary.org/obo/MONDO_0011584 -->
            
            The owl file can be downloaded from https://mondo.monarchinitiative.org/pages/download/
        """

        # Detect file format
        file_format = re.sub(r".*\.", "", mondo_file)

        if file_format != "owl":
            sys.exit("ERROR: Invalid Mondo file format. Accepted format is: owl")

        # Fetch Mondo version from file
        mondo_version = fetch_mondo_version(mondo_file)
        if mondo_version is None:
            sys.exit(f"ERROR: could not detect Mondo version from input file '{mondo_file}'")

        mondo_version_obj = datetime.datetime.strptime(mondo_version, "%Y-%m-%d")

        # Compare G2P OMIM last update with Ensembl (OMIM) data version
        if "Mondo" in g2p_meta_info:
            g2p_mondo_version_obj = datetime.datetime.strptime(g2p_meta_info["Mondo"]["data_version"], "%Y-%m-%d")
            if g2p_mondo_version_obj >= mondo_version_obj:
                sys.exit("Mondo gene-disease data is the latest. Skipping update.")

        # Split the input file - one file per mondo ID
        # This is necessary to correctly map the mondo ID to the hgnc ID
        split_mondo_file(mondo_file)

        # Retrive Mondo gene-disease from Mondo owl file
        print("Getting Mondo gene-disease associations...")
        mondo_gene_diseases, mondo_all_diseases = get_mondo_gene_diseases(mondo_file)
        print("Getting Mondo gene-disease associations... done")

        # We store all Mondo IDs in the table 'disease_external'
        # These IDs are not linked to any gene or other data
        print("> Populating table 'disease_external' with all Mondo IDs...")
        populate_mondo_all_diseases(g2p_db_host, g2p_db_port, g2p_db_name, g2p_user, g2p_password, mondo_all_diseases)
        print("> Populating table 'disease_external' with all Mondo IDs... done")

        # Run the import
        print("> Running Mondo import")
        # Format data for the import
        new_mondo_data = {}
        for mondo, data in mondo_gene_diseases.items():
            if data["hgnc_id"] not in new_mondo_data:
                new_mondo_data[data["hgnc_id"]] = [{ "mondo_id": mondo, "disease": data["disease"] }]
            else:
                new_mondo_data[data["hgnc_id"]].append({ "mondo_id": mondo, "disease": data["disease"] })

        # Insert Mondo gene-disease data into G2P db
        print("Inserting Mondo gene-disease associations into G2P...")
        insert_mondo_gene_diseases(g2p_db_host, g2p_db_port, g2p_db_name, g2p_user, g2p_password, new_mondo_data, mondo_version)
        print("Inserting Mondo gene-disease associations into G2P... done")

if __name__ == '__main__':
    main()
