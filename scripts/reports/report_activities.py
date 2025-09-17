#!/usr/bin/env python3

import sys
import argparse
import configparser
import csv
from datetime import datetime, timedelta
import MySQLdb
import os.path
import requests


"""
    Description: Script to generate a report of changed data in G2P

    Params:
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
    api_url: str, cookies: requests.cookies.RequestsCookieJar
) -> list:
    """
    Fetch the activity logs of the last 7 days.

    Args:
        api_url (str): G2P API URL

    Returns:
        list: list of logs
    """
    url = f"{api_url.rstrip('/')}/activity_logs/?date_cutoff="
    activity_list = []

    # Get the date from the last week
    seven_days_ago = datetime.now().date() - timedelta(days=7)
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
    api_url: str, stable_id:str, cookies: requests.cookies.RequestsCookieJar
) -> list:
    """
    Fetch the activity logs for the specific record.

    Args:
        api_url (str): G2P API URL

    Returns:
        list: list of logs
    """
    url = f"{api_url.rstrip('/')}/activity_logs/?stable_id="
    activity_list = []

    # Get the date from the last week
    seven_days_ago = datetime.now().date() - timedelta(days=7)
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


def format_activity_logs(activity_logs: list) -> dict:
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


def generate_report(output_dir: str, activity_logs_by_record: list):
    current_date = datetime.now().date()
    report_file = "updates_" + str(current_date) + ".txt"

    with (open(report_file, "w") as wr):
        wr.write("Updates from last 7 days\n")
        for g2p_id, list_logs in activity_logs_by_record.items():
            wr.write(f"\n### {g2p_id} ###\n")

            for log in list_logs:
                report_row = None

                # Records logs
                if log["data_type"] == "record":
                    if log["change_type"] == "updated":
                        if log["is_deleted"] == 1:
                            report_row = f"On {log['date']} {log['user']} deleted record {log['g2p_id']}: {log['disease']}; {log['genotype']}; {log['mechanism']}; {log['confidence']}\n"
                        else:
                            # TODO: check what actually was updated
                            # Get the current record data and compare with log
                            report_row = f"(2) On {log['date']} {log['user']} updated record {log['g2p_id']}\n"
                    else:
                        report_row = f"On {log['date']} {log['user']} {log['change_type']} record {log['g2p_id']}: {log['disease']}; {log['genotype']}; {log['mechanism']}; {log['confidence']}\n"
                # Data linked to the record logs
                else:
                    if log["change_type"] == "updated" and "is_deleted" in log and log["is_deleted"] == 1:
                        # TODO: print more details
                        report_row = f"On {log['date']} {log['user']} deleted a {log['data_type']} for record {log['g2p_id']}\n"
                    else:
                        report_row = f"On {log['date']} {log['user']} {log['change_type']} a {log['data_type']} for record {log['g2p_id']}\n"

                wr.write(report_row)


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

    print("Logging in...")
    cookies = login(api_username, api_password, api_url)

    activity_logs = get_activity_logs(api_url, cookies)

    activity_logs_by_record = format_activity_logs(activity_logs)

    generate_report(output_dir, activity_logs_by_record)

    print("Logging out...")
    logout(api_url, cookies)


if __name__ == "__main__":
    main()
