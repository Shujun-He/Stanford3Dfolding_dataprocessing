import requests
import json
from multiprocessing import Pool, cpu_count
import os

# Directory to save the PDB, CIF, and FASTA files
SAVE_DIR = "rna_structures"
os.makedirs(SAVE_DIR, exist_ok=True)

# Base URLs for RCSB PDB APIs
SEARCH_API = "https://search.rcsb.org/rcsbsearch/v2/query"
DOWNLOAD_PDB = "https://files.rcsb.org/download/{pdb_id}.pdb"
DOWNLOAD_CIF = "https://files.rcsb.org/download/{pdb_id}.cif"
DOWNLOAD_FASTA = "https://www.rcsb.org/fasta/entry/{pdb_id}"
SUMMARY_API = "https://data.rcsb.org/rest/v1/core/structure/{pdb_id}"

# Query to search for RNA structures in PDB
def search_rna_structures():
    query = {
        "query": {
            "type": "terminal",
            "service": "full_text",
            "parameters": {
                "value": "RNA"
            }
        },
        "return_type": "entry",
        "request_options": {
            "return_all_hits": True,
            "results_content_type": ["experimental"],
            "sort": [{"sort_by": "score", "direction": "desc"}]
        }
    }
    response = requests.post(SEARCH_API, json=query)
    if response.status_code == 200:
        return [entry['identifier'] for entry in response.json()['result_set']]
    else:
        print("Error in search API:", response.text)
        return []

# Fetch resolution and save PDB/CIF/FASTA files
def fetch_structure_data(pdb_id):
    print(f"Processing PDB ID: {pdb_id}")
    structure_info = {"PDB_ID": pdb_id, "Resolution": "N/A", "PDB_File": None, "CIF_File": None, "FASTA_File": None}

    # Fetch resolution and metadata
    response = requests.get(SUMMARY_API.format(pdb_id=pdb_id))
    if response.status_code == 200:
        data = response.json()
        structure_info["Resolution"] = data.get("rcsb_entry_info", {}).get("resolution_combined", ["N/A"])[0]

    # Download PDB file
    pdb_file_path = os.path.join(SAVE_DIR, f"{pdb_id}.pdb")
    if download_and_save(DOWNLOAD_PDB.format(pdb_id=pdb_id), pdb_file_path):
        structure_info["PDB_File"] = pdb_file_path

    # Download CIF file
    cif_file_path = os.path.join(SAVE_DIR, f"{pdb_id}.cif")
    if download_and_save(DOWNLOAD_CIF.format(pdb_id=pdb_id), cif_file_path):
        structure_info["CIF_File"] = cif_file_path

    # Download FASTA file
    fasta_file_path = os.path.join(SAVE_DIR, f"{pdb_id}.fasta")
    if download_and_save(DOWNLOAD_FASTA.format(pdb_id=pdb_id), fasta_file_path):
        structure_info["FASTA_File"] = fasta_file_path

    return structure_info

# Helper function to download and save files
def download_and_save(url, file_path):
    try:
        response = requests.get(url)
        if response.status_code == 200:
            with open(file_path, "wb") as file:
                file.write(response.content)
            return True
        else:
            print(f"Failed to download {url}")
            return False
    except Exception as e:
        print(f"Error downloading {url}: {e}")
        return False

# Main execution
def main():
    print("Searching for RNA structures...")
    rna_pdb_ids = search_rna_structures()
    print(f"Found {len(rna_pdb_ids)} RNA structures.")
    # exit()
    print("Fetching structure data using multiprocessing...")
    with Pool(cpu_count()) as pool:
        structure_data = pool.map(fetch_structure_data, rna_pdb_ids)

    # Save metadata to JSON file
    metadata_file = os.path.join(SAVE_DIR, "rna_structures_metadata.json")
    with open(metadata_file, "w") as json_file:
        json.dump(structure_data, json_file, indent=4)
    print(f"Saved structure metadata to {metadata_file}")

if __name__ == "__main__":
    main()