#!/usr/bin/env python3

import argparse
import os.path
import sys
from datetime import datetime, timedelta

import requests
from fpdf import FPDF
from fpdf.enums import XPos, YPos


"""
    Description: Script to generate a report of data updates from the last 7 days.
                 It writes two report files:
                    - report_all_updates_YYYY-MM-DD.pdf
                    - report_updates_YYYY-MM-DD.pdf

    Options:
            --api-url      : G2P API URL (mandatory)
            --api-username : Username to connect to the G2P API (mandatory)
            --api-password : Password to connect to the G2P API (mandatory)
            --output-dir   : Path to the output directory where the report is going to be saved (mandatory)
    """


def login(
    api_username: str, api_password: str, api_url: str
) -> requests.cookies.RequestsCookieJar:
    """Login into G2P API"""
    login_url = f"{api_url.rstrip('/')}/login/"

    response = requests.post(
        login_url, json={"username": api_username, "password": api_password}
    )

    if response.status_code != 200:
        sys.exit("Login failed. Check your credentials and API URL.")

    return response.cookies


def logout(api_url: str, cookies: requests.cookies.RequestsCookieJar) -> None:
    """Logout of the API"""
    logout_url = f"{api_url.rstrip('/')}/logout/"

    response = requests.post(logout_url, cookies=cookies)

    if response.status_code != 204:
        sys.exit("Logout failed. Check your credentials and API URL.")


def get_activity_logs(
    api_url: str, seven_days_ago: datetime, cookies: requests.cookies.RequestsCookieJar
) -> list:
    """
    Fetch the activity logs of the last 7 days.

    Args:
        api_url (str): G2P API URL
        seven_days_ago (datetime): date from 7 days ago
        cookies (requests.cookies.RequestsCookieJar): cookies containing the login token

    Returns:
        list: list of logs
    """
    url = f"{api_url.rstrip('/')}/activity_logs/?date_cutoff="
    activity_list = []

    url_last_week = url + str(seven_days_ago)

    while url_last_week:
        try:
            response = requests.get(url_last_week, cookies=cookies)
        except Exception as e:
            sys.exit("Error while fetching activity logs:", str(e))
        else:
            if response.status_code == 200:
                result = response.json()
                activity_list += result["results"]
                url_last_week = result["next"]
            else:
                sys.exit(
                    f"Failed to fetch activity logs from API. Status code: {response.status_code}"
                )

    return activity_list


def get_record_activity_logs(
    api_url: str, stable_id: str, cookies: requests.cookies.RequestsCookieJar
) -> list:
    """
    Fetch the activity logs for the specific record.

    Args:
        api_url (str): G2P API URL
        stable_id (str): G2P stable ID of the record
        cookies (requests.cookies.RequestsCookieJar): cookies containing the login token

    Returns:
        list: list of logs
    """
    url_record_activity = f"{api_url.rstrip('/')}/activity_logs/?stable_id={stable_id}"
    activity_list = []

    while url_record_activity:
        try:
            response = requests.get(url_record_activity, cookies=cookies)
        except Exception as e:
            sys.exit("Error while fetching activity logs:", str(e))
        else:
            if response.status_code == 200:
                result = response.json()
                for log in result["results"]:
                    # We are only interested in the log type "record"
                    if log["data_type"] == "record":
                        activity_list.append(log)
                url_record_activity = result["next"]
            else:
                sys.exit(
                    f"Failed to fetch activity logs from API for record {stable_id}. Status code: {response.status_code}"
                )

    return activity_list


def get_record_update(activity_list: list, log_date: str) -> list:
    """
    Get the specific data update from the record activity logs.

    Args:
        activity_list (list): list of activity logs
        log_date (str): date of the activity log

    Returns:
        str: specific update done on the date defined in log_date
    """
    update = ""

    for idx, log in enumerate(activity_list):
        if log["date"] == log_date:
            current_data = log
            # Check if there is a previous activity log
            if idx + 1 < len(activity_list):
                previous_data = activity_list[idx + 1]

                # Check what is different
                if current_data["confidence"] != previous_data["confidence"]:
                    update += f" from confidence '{previous_data['confidence']}' to '{current_data['confidence']}';"
                if current_data["mechanism"] != previous_data["mechanism"]:
                    update += f" from mechanism '{previous_data['mechanism']}' to '{current_data['mechanism']}';"
                if (
                    current_data["mechanism_support"]
                    != previous_data["mechanism_support"]
                ):
                    update += f" from mechanism support '{previous_data['mechanism_support']}' to '{current_data['mechanism_support']}';"
                if current_data["disease"] != previous_data["disease"]:
                    update += f" from disease '{previous_data['disease']}' to '{current_data['disease']}';"
                if current_data["is_reviewed"] != previous_data["is_reviewed"]:
                    update += f" from is_reviewed '{previous_data['is_reviewed']}' to '{current_data['is_reviewed']}';"
                if current_data["is_deleted"] != previous_data["is_deleted"]:
                    update += f" from is_deleted '{previous_data['is_deleted']}' to '{current_data['is_deleted']}';"

    return update


def format_activity_logs(activity_logs: list) -> dict:
    """
    Returns the activity logs by the G2P stable ID.

    Args:
        activity_list (list): list of activity logs

    Returns:
        dict: key is the stable ID and the value is the list of activity logs
    """
    activity_by_record = {}

    for log in activity_logs:
        if "g2p_id" not in log:
            continue
        elif log["g2p_id"] in activity_by_record:
            activity_by_record[log["g2p_id"]].append(log)
        else:
            activity_by_record[log["g2p_id"]] = []
            activity_by_record[log["g2p_id"]].append(log)

    return activity_by_record


def generate_report(
    output_dir: str,
    activity_logs_by_record: list,
    seven_days_ago: datetime,
    api_url: str,
    cookies: requests.cookies.RequestsCookieJar,
) -> None:
    """
    Generate the final reports of all the data updates done in the last 7 days.

    Args:
        output_dir (str): directory where the reports are going to be saved
        activity_logs_by_record (dict): list of activity logs for each G2P stable ID
        seven_days_ago (datetime): date
        api_url (str): G2P API URL
        cookies (requests.cookies.RequestsCookieJar): cookies containing the login token
    """
    current_date = datetime.now().date()
    full_report_file = os.path.join(
        output_dir, "report_all_updates_" + str(current_date) + ".pdf"
    )
    minimal_report_file = os.path.join(
        output_dir, "report_updates_" + str(current_date) + ".pdf"
    )

    # Create the PDF with all the data updates - full report
    pdf = FPDF()
    pdf.add_page()
    pdf.set_font("helvetica", style="B", size=16)
    pdf.multi_cell(0, 7, f"Data updates from {seven_days_ago} to {current_date}\n")
    pdf.set_font("helvetica", size=12)

    # Create the PDF with only the record data updates - minimal report
    pdf_minimal = FPDF()
    pdf_minimal.add_page()
    pdf_minimal.set_font("helvetica", style="B", size=16)
    pdf_minimal.multi_cell(
        0, 7, f"Data updates from {seven_days_ago} to {current_date}\n"
    )
    pdf_minimal.set_font("helvetica", size=12)

    for g2p_id, list_logs in activity_logs_by_record.items():
        pdf.set_text_color(0, 0, 255)
        pdf.multi_cell(
            0,
            5,
            f"\n### {g2p_id} ###\n",
            new_x=XPos.LMARGIN,
            new_y=YPos.NEXT,
            link="https://www.ebi.ac.uk/gene2phenotype/lgd/" + g2p_id,
        )
        pdf.set_text_color(0, 0, 0)
        # Get all the activity logs for the record
        record_activity_logs = get_record_activity_logs(api_url, g2p_id, cookies)
        list_report_minimal = []

        for log in list_logs:
            report_row = None
            report_minimal_row = None

            # Records logs
            if log["data_type"] == "record":
                if log["change_type"] == "updated":
                    if log["is_deleted"] == 1:
                        report_row = f"On {log['date']} {log['user']} deleted record {log['g2p_id']}: {log['disease']}; {log['genotype']}; {log['mechanism']}; {log['confidence']}\n"
                        # Print to the minimal report
                        # pdf_minimal.multi_cell(0, 5, report_row)
                        report_minimal_row = report_row
                    else:
                        # Get the current record data and compare with log
                        record_updates = get_record_update(
                            record_activity_logs, log["date"]
                        )
                        if record_updates != "":
                            report_row = f"On {log['date']} {log['user']} updated record {log['g2p_id']}:{record_updates}\n"
                            # Print to the minimal report
                            # pdf_minimal.multi_cell(0, 5, report_row)
                            report_minimal_row = report_row
                else:
                    report_row = f"On {log['date']} {log['user']} {log['change_type']} record {log['g2p_id']}: {log['disease']}; {log['genotype']}; {log['mechanism']}; {log['confidence']}\n"
                    # Print to the minimal report
                    # pdf_minimal.multi_cell(0, 5, report_row)
                    report_minimal_row = report_row
            # Data linked to the record logs
            else:
                if (
                    log["change_type"] == "updated"
                    and "is_deleted" in log
                    and log["is_deleted"] == 1
                ):
                    report_row = f"On {log['date']} {log['user']} deleted a {log['data_type']} for record {log['g2p_id']}\n"
                elif log["change_type"] == "created":
                    report_row = (
                        f"On {log['date']} {log['user']} created {log['data_type']}: "
                    )
                    if log["data_type"] == "panel":
                        report_row += f"{log['panel_name']}\n"
                        # Print to the minimal report
                        report_minimal_row = f"On {log['date']} {log['user']} created {log['data_type']}: {log['panel_name']} for record {log['g2p_id']}\n"
                        # pdf_minimal.multi_cell(0, 5, f"On {log['date']} {log['user']} created {log['data_type']}: {log['panel_name']} for record {log['g2p_id']}\n")
                    if log["data_type"] == "publication":
                        report_row += f"PMID {log['publication_pmid']}\n"
                    if log["data_type"] == "phenotype":
                        report_row += (
                            f"{log['phenotype']} for PMID {log['publication_pmid']}\n"
                        )
                    if log["data_type"] == "phenotype_summary":
                        report_row += f"'{log['summary']}'\n"
                    if log["data_type"] == "variant_consequence":
                        report_row += f"'{log['variant_consequence']}'\n"
                    if log["data_type"] == "variant_type":
                        report_row += f"'{log['variant_type']}' for PMID {log['publication_pmid']}\n"
                    if log["data_type"] == "record_comment":
                        is_public = (
                            "public comment" if log["is_public"] else "private comment"
                        )
                        report_row += f"{is_public}\n"
                    if log["data_type"] == "mechanism_synopsis":
                        report_row += (
                            f"'{log['synopsis']}' with support '{log['support']}'\n"
                        )
                    if log["data_type"] == "mechanism_evidence":
                        report_row += f"evidence type '{log['evidence_type']}', evidence value '{log['evidence']}' for PMID {log['publication_pmid']}\n"
                    if log["data_type"] == "cross_cutting_modifier":
                        report_row += f"'{log['ccm']}'\n"
                else:
                    report_row = f"On {log['date']} {log['user']} {log['change_type']} a {log['data_type']} for record {log['g2p_id']}\n"

            # Write to the full report
            if report_row:
                pdf.multi_cell(0, 5, report_row)

            # Save which rows are going to be reported in the minimal report
            if report_minimal_row:
                list_report_minimal.append(report_minimal_row)

        # For each G2P ID write to the minimal report
        if list_report_minimal:
            pdf_minimal.set_text_color(0, 0, 255)
            pdf_minimal.multi_cell(
                0,
                5,
                f"\n### {g2p_id} ###\n",
                new_x=XPos.LMARGIN,
                new_y=YPos.NEXT,
                link="https://www.ebi.ac.uk/gene2phenotype/lgd/" + g2p_id,
            )
            pdf_minimal.set_text_color(0, 0, 0)
            for row in list_report_minimal:
                pdf_minimal.multi_cell(0, 5, row)

    pdf.output(full_report_file)
    pdf_minimal.output(minimal_report_file)


def main():
    parser = argparse.ArgumentParser(
        description="Generate a G2P XML file for the EBI search engine"
    )
    parser.add_argument("--api-url", required=True, help="G2P API URL")
    parser.add_argument(
        "--api-username", required=True, help="Username to connect to the G2P API"
    )
    parser.add_argument(
        "--api-password", required=True, help="Password to connect to the G2P API"
    )
    parser.add_argument(
        "--output-dir", required=True, help="Path to the output directory"
    )
    args = parser.parse_args()

    api_url = args.api_url
    api_username = args.api_username
    api_password = args.api_password
    output_dir = args.output_dir

    if not os.path.exists(output_dir):
        sys.exit(f"Invalid output directory '{output_dir}'")

    # Get the date from the last week
    seven_days_ago = datetime.now().date() - timedelta(days=7)

    print("Logging in...")
    cookies = login(api_username, api_password, api_url)

    print("Getting the activity logs from the API...")
    activity_logs = get_activity_logs(api_url, seven_days_ago, cookies)
    activity_logs_by_record = format_activity_logs(activity_logs)
    print("Getting the activity logs from the API... done")

    print("Generating reports...")
    generate_report(
        output_dir, activity_logs_by_record, seven_days_ago, api_url, cookies
    )
    print("Generating reports... done")

    print("Logging out...")
    logout(api_url, cookies)


if __name__ == "__main__":
    main()
