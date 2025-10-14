from pathlib import Path
import requests
import json
import os
from functools import partial
from retryhttp import retry
from utils import parallel_process

# Base URLs for RCSB PDB APIs
SEARCH_API = "https://search.rcsb.org/rcsbsearch/v2/query"
DOWNLOAD_PDB = "https://files.rcsb.org/download/{pdb_id}.pdb"
DOWNLOAD_CIF = "https://files.rcsb.org/download/{pdb_id}.cif"
DOWNLOAD_FASTA = "https://www.rcsb.org/fasta/entry/{pdb_id}"
SUMMARY_API = "https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"


SEARCH_QUERIES = {
    "full_text": {
        "type": "terminal",
        "service": "full_text",
        "parameters": {"value": "__UNDEFINED__"},
    },
    "polymer_type": {
        "type": "terminal",
        "service": "text",
        "parameters": {
            "operator": "exact_match",
            "value": "__UNDEFINED__",
            "attribute": "entity_poly.rcsb_entity_polymer_type",
        },
    },
}


# Query to search for RNA structures in PDB, limit dates and fetch only structures with RNA entities
def search_rna_structures(
    start_date="now-1w", stop_date="now", search_type="polymer_type", search_term="RNA"
):
    search_query = SEARCH_QUERIES.get(search_type)
    if search_query is None:
        raise ValueError(f"Unsupported search_type: {search_type}")
    search_query["parameters"]["value"] = search_term

    query = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_accession_info.initial_release_date",
                        "operator": "greater",
                        "value": f"{start_date}",
                    },
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_accession_info.initial_release_date",
                        "operator": "less_or_equal",
                        "value": f"{stop_date}",
                    },
                },
                SEARCH_QUERIES[search_type],
            ],
        },
        "return_type": "entry",
        "request_options": {
            "return_all_hits": True,
            "results_content_type": ["experimental"],
            "sort": [{"sort_by": "score", "direction": "desc"}],
        },
    }
    response = requests.post(SEARCH_API, json=query)
    if response.status_code == 200:
        return [entry["identifier"] for entry in response.json()["result_set"]]
    else:
        print("Error in search API:", response.text)
        return []


# Extra metadata fields to extract from the summary API, key is the desired field name, value is dot notation key to extract
extra_metadata_fields = {
    "Release": "rcsb_accession_info.initial_release_date",
    "Resolution": "rcsb_entry_info.resolution_combined.0",
    "Method": "rcsb_entry_info.experimental_method",
    "Title": "struct.title",
    "Description": "struct.title",
    "Keywords": "struct_keywords.pdbx_keywords",
}


# TODO: Use jsonpath for more complex queries if needed
def get_dot_notation_item(data, path):
    keys = path.split(".")
    value = data
    current_path = ""
    for key in keys:
        # Check if map or list, and access accordingly
        if isinstance(value, dict):

            def accessor(v, k):
                return v.get(k, None)
        elif isinstance(value, list):
            # Convert key to int
            def accessor(v, k):
                return v[int(k)]
        else:
            raise ValueError(
                f"Data at path {current_path} is neither dict nor list. Unable to access key {key} "
            )
        value = accessor(value, key)
        if value is None:
            raise ValueError(f"Key {key} not found in data[{current_path}]")
        current_path += f".{key}" if current_path else key
    return value


# Function to extract specified metadata fields from the JSON response
def get_metadata_fields(metadata_file, metadata_fields):
    metadata = {}
    with open(metadata_file, "r") as f:
        data = json.load(f)
        for field, path in metadata_fields.items():
            try:
                value = get_dot_notation_item(data, path)
                metadata[field] = value
            except Exception as e:
                print(f"Error reading {path} from {metadata_file}: {e}")
                metadata[field] = None
    return metadata


# Fetch metadata and save PDB/CIF/FASTA files
def fetch_structure_data(pdb_id, output_dir):
    print(f"Processing PDB ID: {pdb_id}")
    structure_info = {
        "PDB_ID": pdb_id,
        "PDB_File": None,
        "CIF_File": None,
        "FASTA_File": None,
        "Metadata_File": None,
    }

    # Fetch metadata and store it in json, do not parse
    metadata_path = os.path.join(output_dir, f"{pdb_id}.metadata.json")
    if download_and_save(SUMMARY_API.format(pdb_id=pdb_id), metadata_path):
        structure_info["Metadata_File"] = metadata_path
    # Extract specified metadata fields
    metadata = get_metadata_fields(metadata_path, extra_metadata_fields)
    structure_info.update(metadata)

    # Download PDB file
    pdb_file_path = os.path.join(output_dir, f"{pdb_id}.pdb")
    if download_and_save(DOWNLOAD_PDB.format(pdb_id=pdb_id), pdb_file_path):
        structure_info["PDB_File"] = pdb_file_path

    # Download CIF file
    cif_file_path = os.path.join(output_dir, f"{pdb_id}.cif")
    if download_and_save(DOWNLOAD_CIF.format(pdb_id=pdb_id), cif_file_path):
        structure_info["CIF_File"] = cif_file_path

    # Download FASTA file
    fasta_file_path = os.path.join(output_dir, f"{pdb_id}.fasta")
    if download_and_save(DOWNLOAD_FASTA.format(pdb_id=pdb_id), fasta_file_path):
        structure_info["FASTA_File"] = fasta_file_path

    return structure_info

def is_retryable_http_error(result):
    if isinstance(result, requests.Response):
        return result.status_code == 429  # Retry only on 429
    return False

# Helper function to download and save files,
@retry
def download_and_save(url, file_path, force=False):
    """Download a file from a URL and save it locally. Skip if file exists unless force is True.

    Retry if the download fails due to 429 Too Many Requests or other transient errors."""
    if Path(file_path).exists() and not force:
        print(f"File '{file_path}' already exists. Skipping download.")
        return True
    try:
        response = requests.get(url)
        response.raise_for_status()
        with open(file_path, "wb") as file:
            file.write(response.content)
        return True
    except requests.exceptions.RequestException as err:
        print(f"Failed to download '{url}': {err}")
        return False


# Main execution
def main():
    import argparse

    parser = argparse.ArgumentParser(description="Fetch RNA structures from PDB")
    parser.add_argument(
        "--start_date",
        type=str,
        default="now-1w",
        help="Start date for search (e.g., 'now-1w')",
    )
    parser.add_argument(
        "--stop_date",
        type=str,
        default="now",
        help="Stop date for search (e.g., 'now')",
    )
    parser.add_argument(
        "--search_type",
        type=str,
        choices=SEARCH_QUERIES.keys(),
        default="polymer_type",
        help="Type of search query",
    )
    parser.add_argument(
        "--search_term", type=str, default="RNA", help="Term to search for"
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        default="rna_structures_all",
        help="Directory to save downloaded files",
    )
    args = parser.parse_args()

    print(
        f"Searching for RNA structures released from {args.start_date} to {args.stop_date}..."
    )
    print(f"Using search type: {args.search_type}")
    rna_pdb_ids = search_rna_structures(
        start_date=args.start_date,
        stop_date=args.stop_date,
        search_type=args.search_type,
        search_term=args.search_term,
    )
    print(f"Found {len(rna_pdb_ids)} RNA structures.")
    if len(rna_pdb_ids) > 0:
        os.makedirs(args.output_dir, exist_ok=True)

    print("Fetching structure data using multiprocessing...")
    structure_data = parallel_process(
        partial(fetch_structure_data, output_dir=args.output_dir),
        rna_pdb_ids,
        num_processes=4,  # Throttle to limit number of requests
        desc="Fetching structure data",
    )

    # Save metadata to JSON file
    metadata_file = os.path.join(args.output_dir, "rna_structures_metadata.json")
    with open(metadata_file, "w") as json_file:
        json.dump(structure_data, json_file, indent=4)
    print(f"Saved structure metadata to {metadata_file}")


if __name__ == "__main__":
    main()
