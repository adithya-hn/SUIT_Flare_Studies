#####################################################
# @Author: SoLEXSPOC
# @Date:   2025-10-25 20:19:09 pm
# @email: sarwade@ursc.gov.in
# @File Name: caldb_utils.py
# @Project: solexs_tools
#
# @Last Modified time: 2026-01-06 05:37:55 pm
#####################################################

import argparse
import importlib.resources
import shutil, os, json, datetime
from pathlib import Path

def get_caldb_base_dir():
    CALDB_BASE_DIR = os.environ.get('SOLEXS_CALDB')
    
    if CALDB_BASE_DIR is None:
        raise EnvironmentError(
            "\n\n*** ERROR: SOLEXS_CALDB environment variable is not set. ***\n"
            "Please set this variable to point to your unzipped CALDB directory.\n"
            "Example: export SOLEXS_CALDB=\"/path/to/my/solexs_caldb\"\n"
            "To extract the default CALDB, run this command:\n"
            "solexs-caldb-extract /path/to/my/solexs_caldb\n"
        )            
    return CALDB_BASE_DIR


def walk_and_copy_resources(source_path, target_path: Path):
    for resource in source_path.iterdir():
        resource_name = resource.name
        
        dest_item = target_path / resource_name
        
        if resource.is_dir():
            dest_item.mkdir(parents=True, exist_ok=True)
            walk_and_copy_resources(resource, dest_item)
        elif resource.is_file():
            try:
                with importlib.resources.as_file(resource) as src:
                    shutil.copyfile(src, dest_item)
                    print(f"Copied: {dest_item.relative_to(target_path.parent)}")
            except Exception as e:
                print(f"Warning: Could not copy {resource_name}. Error: {e}")

def extract_caldb(target_path: Path):
    """
    Extracts the packaged CALDB to a specified target directory.
    """
    
    try:
        source_caldb_root = importlib.resources.files('solexs_tools').joinpath('CALDB')
    except (ImportError, AttributeError):
        print("Error: Could not find package data. Make sure solexs_tools is installed.")
        return

    print(f"Source CALDB found. Extracting to: {target_path}")
    
    try:
        target_path.mkdir(parents=True, exist_ok=True)
        
        walk_and_copy_resources(source_caldb_root, target_path)

        print("\nExtraction complete.")
        print("To use this CALDB, set your environment variable:")
        print(f"\n  export SOLEXS_CALDB=\"{target_path.resolve()}\"\n")
        print("Add this line to your ~/.bashrc or ~/.zshrc file to make it permanent.")
        
    except Exception as e:
        print(f"\nAn error occurred during extraction: {e}")
        print("Please check file permissions and try again.")


def solexs_caldb_extract_cli():
    parser = argparse.ArgumentParser(
        description="Extract the default SoLEXS CALDB to a specified directory."
    )
    parser.add_argument(
        "output_directory", 
        type=str, 
        help="The target directory to copy CALDB files (e.g., /home/user/solexs_caldb)"
    )
    args = parser.parse_args()
    
    target_path = Path(args.output_directory).resolve()
    
    if target_path.exists() and not target_path.is_dir():
        print(f"Error: Target path '{target_path}' exists and is not a directory.")
        return

    if target_path.exists() and any(target_path.iterdir()):
        print(f"Warning: Directory '{target_path}' is not empty.")
        confirm = input("Do you want to proceed and potentially overwrite files? (y/n): ")
        if confirm.lower() != 'y':
            print("Extraction cancelled.")
            return

    try:
        extract_caldb(target_path)
    except Exception as e:
        print(f"An error occurred: {e}")


def get_caldb_file(file_type, filter_sdd='SDD2', obs_date=None):
    
    CALDB_BASE_DIR = get_caldb_base_dir()

    index_path = os.path.join(CALDB_BASE_DIR, "caldb_index.json")
    if not os.path.exists(index_path):
        raise FileNotFoundError(f"CALDB Index missing: {index_path}")

    with open(index_path, 'r') as f:
        index = json.load(f)


    det_entries = index['files'].get(filter_sdd)
    if det_entries is None:
        raise KeyError(f"Detector '{filter_sdd}' not found in CALDB index.")

    entries = det_entries.get(file_type)
    if not entries:
        raise ValueError(f"No valid '{file_type}' files defined for detector '{filter_sdd}'.")

    target_dt = datetime.datetime.utcnow()
    
    if obs_date is not None:
        target_dt = datetime.datetime.strptime(obs_date, '%Y%m%d')

    selected_file = None
    
    for entry in entries:
        start_dt = datetime.datetime.fromisoformat(entry['valid_start'])
        
        if entry['valid_stop']:
            stop_dt = datetime.datetime.fromisoformat(entry['valid_stop'])
        else:
            stop_dt = datetime.datetime.max 

        if start_dt <= target_dt < stop_dt:
            selected_file = entry['filename']
            break
    
    if not selected_file:
        raise ValueError(f"No {file_type} file found for {filter_sdd} at time {target_dt}")

    return os.path.join(CALDB_BASE_DIR, selected_file)