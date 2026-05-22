#####################################################
# @Author: SoLEXSPOC
# @Date:   2024-11-17 11:44:26 pm
# @email: sarwade@ursc.gov.in
# @File Name: time_utils.py
# @Project: solexs_tools
#
# @Last Modified time: 2025-01-29 04:55:06 pm
#####################################################

import argparse
from datetime import datetime, timezone

def unix_time_to_utc(unix_time):
    """
    Convert Unix timestamp to UTC in ISO 8601 format.

    Args:
        unix_time (int): Unix timestamp.

    Returns:
        str: UTC time in ISO 8601 format.
    """
    try:
        dt = datetime.fromtimestamp(unix_time, timezone.utc)
        return dt.isoformat()
    except Exception as e:
        raise ValueError(f"Invalid Unix time: {unix_time}. Error: {e}")


def utc_to_unix_time(utc_time_str):
    """
    Convert UTC time in ISO 8601 format to Unix timestamp.

    Args:
        utc_time_str (str): UTC time in ISO 8601 format.

    Returns:
        int: Unix timestamp.
    """
    try:
        dt = datetime.fromisoformat(utc_time_str + '+00:00' )
        return int(dt.timestamp())
    except Exception as e:
        raise ValueError(f"Invalid UTC time format: {utc_time_str}. Error: {e}")


def solexs_time2utc_cli():
    parser = argparse.ArgumentParser(
        description="Convert Unix timestamp to UTC in ISO 8601 format."
    )
    parser.add_argument("unix_time", type=int, help="Unix timestamp (e.g., 1633046400)")
    args = parser.parse_args()

    try:
        utc_time = unix_time_to_utc(args.unix_time)
        print(f"UTC Time: {utc_time}")
    except ValueError as e:
        print(e)

def solexs_utc2time_cli():
    parser = argparse.ArgumentParser(
        description="Convert UTC in ISO 8601 format to Unix timestamp."
    )
    parser.add_argument(
        "utc_time", type=str, help="UTC time in ISO 8601 format (e.g., 2021-10-01T00:00:00)"
    )
    args = parser.parse_args()

    try:
        unix_time = utc_to_unix_time(args.utc_time)
        print(f"Unix Timestamp: {unix_time}")
    except ValueError as e:
        print(e)