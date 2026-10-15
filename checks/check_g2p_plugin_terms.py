#!/usr/bin/env python3

import argparse
import configparser
from email.message import EmailMessage
import json
from pathlib import Path
import re
import smtplib
import sys
from typing import Optional, Sequence, TypedDict
from urllib.error import URLError
from urllib.request import urlopen

"""
Description:
    Check the G2P plugin's supported terms against the G2P attributes API.
    The plugin source is downloaded from Ensembl/VEP_plugins on GitHub (main).
    Legacy plugin terms absent from the G2P API are allowed.
    The unsupported confidence levels disputed and refuted are ignored.
    This script only flags missing terms, it never changes the plugin.

How to run the script:
    python scripts/check_g2p_terms.py --email-config /path/to/config.ini

Output:
    0 = covered
    1 = missing terms
    2 = check, configuration or email delivery failed
"""


API_URL = "https://www.ebi.ac.uk/gene2phenotype/api/attribs/"
PLUGIN_URL = "https://raw.githubusercontent.com/Ensembl/VEP_plugins/main/G2P.pm"
TERM_TABLES = {
    "genotype": "allelic_requirements",
    "confidence_category": "supported_confidence_levels",
}
IGNORED_CONFIDENCE_LEVELS = {"disputed", "refuted"}


class EmailConfig(TypedDict):
    sender: str
    host: str
    recipients: list[str]
    enabled: bool
    port: int


def plugin_terms(source: str, table: str) -> set[str]:
    """Read literal keys without loading Perl or its Ensembl dependencies.

    Fail if the table format changes instead of silently reporting coverage.
    """
    blocks = re.findall(
        r"^my\s+\$" + re.escape(table) + r"\s*=\s*\{\s*\n(.*?)^\};",
        source,
        re.MULTILINE | re.DOTALL,
    )
    if len(blocks) != 1:
        raise ValueError(f"Expected one literal ${table} table in the plugin")
    terms = set()
    for line in blocks[0].splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        match = re.match(r"\s*'([^'\\]+)'\s*=>", line)
        if not match:
            raise ValueError(f"Cannot parse ${table} entry: {line.strip()}")
        if match[1] in terms:
            raise ValueError(f"Duplicate ${table} term: {match[1]}")
        terms.add(match[1])
    if not terms:
        raise ValueError(f"Empty ${table} table in the plugin")
    return terms


def missing_terms(source: str, attributes: object) -> dict[str, list[str]]:
    if not isinstance(attributes, dict):
        raise ValueError("API response must be a JSON object")
    missing = {}
    for field, table in TERM_TABLES.items():
        terms = attributes.get(field)
        if (
            not isinstance(terms, list)
            or not terms
            or any(not isinstance(term, str) or not term.strip() for term in terms)
        ):
            raise ValueError(f"API field {field!r} must be a nonempty list of terms")
        required_terms = set(terms)
        if field == "confidence_category":
            required_terms -= IGNORED_CONFIDENCE_LEVELS
        missing[field] = sorted(required_terms - plugin_terms(source, table))
    return missing


def positive_timeout(value: str) -> int:
    timeout = int(value)
    if timeout <= 0:
        raise argparse.ArgumentTypeError("timeout must be a positive integer")
    return timeout


def read_email_config(path: Path) -> EmailConfig:
    """Read the [email] section of the EBI SMTP configuration."""
    parser = configparser.ConfigParser(interpolation=None)
    parser.read_string(path.read_text(encoding="utf-8"))
    if not parser.has_section("email"):
        raise ValueError("Missing [email] section in email configuration")
    config = parser["email"]
    for key in ("from", "host", "mailing_list", "send_to_mailing_list", "port"):
        if not config.get(key, "").strip():
            raise ValueError(f"Missing email configuration value: {key}")
    enabled = config.getboolean("send_to_mailing_list")
    port = config.getint("port")
    if not 1 <= port <= 65535:
        raise ValueError("SMTP port must be between 1 and 65535")
    recipients = [
        address.strip()
        for address in config["mailing_list"].split(",")
        if address.strip()
    ]
    if not recipients:
        raise ValueError("Email mailing_list must contain a recipient")
    return {
        "sender": config["from"].strip(),
        "host": config["host"].strip(),
        "recipients": recipients,
        "enabled": enabled,
        "port": port,
    }


def send_missing_terms_email(
    config: EmailConfig, report: str, plugin: str, url: str, timeout: int
) -> None:
    """Use the EBI SMTP relay without authentication, as configured."""
    message = EmailMessage()
    message["From"] = config["sender"]
    message["To"] = ", ".join(config["recipients"])
    message["Subject"] = "G2P plugin: missing supported terms"
    message.set_content(f"Plugin: {plugin}\nG2P API: {url}\n\n{report}\n")
    with smtplib.SMTP(config["host"], config["port"], timeout=timeout) as smtp:
        refused = smtp.send_message(message)
        if refused:
            raise smtplib.SMTPRecipientsRefused(refused)


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--url", default=API_URL, help="G2P API attribs endpoint")
    parser.add_argument(
        "--timeout", type=positive_timeout, default=30, help="HTTP timeout in seconds"
    )
    parser.add_argument(
        "--quiet", action="store_true", help="print only missing terms or errors"
    )
    parser.add_argument(
        "--email-config",
        type=Path,
        help="path to the EBI SMTP INI configuration file with an [email] section",
    )
    args = parser.parse_args(argv)
    try:
        email_config = (
            read_email_config(args.email_config) if args.email_config else None
        )
        with urlopen(PLUGIN_URL, timeout=args.timeout) as response:
            source = response.read().decode("utf-8")
        with urlopen(args.url, timeout=args.timeout) as response:
            attributes = json.load(response)
        missing = missing_terms(source, attributes)
    except (OSError, URLError, ValueError, configparser.Error) as error:
        print(f"G2P term check failed: {error}", file=sys.stderr)
        return 2

    report_lines = []
    for field, terms in missing.items():
        if terms:
            report_lines.append(
                f"Missing {field} in ${TERM_TABLES[field]}: {', '.join(terms)}"
            )
    if any(missing.values()):
        report_lines.append(
            "Review these terms and update the G2P plugin where appropriate."
        )
        report = "\n".join(report_lines)
        print(report)
        if email_config and email_config["enabled"]:
            try:
                send_missing_terms_email(
                    email_config, report, PLUGIN_URL, args.url, args.timeout
                )
            except (OSError, smtplib.SMTPException, ValueError) as error:
                print(f"G2P email notification failed: {error}", file=sys.stderr)
                return 2
        return 1
    if not args.quiet:
        print(
            "OK: G2P plugin covers all required G2P genotype and confidence terms (disputed and refuted excluded)."
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
