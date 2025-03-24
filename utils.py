import requests

# Function to fetch the publication date for a PDB ID
def get_pdb_publication_date(pdb_id):
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id.lower()}"
    try:
        response = requests.get(url)
        response.raise_for_status()
        data = response.json()
        # Extract the deposition date
        return data.get('rcsb_accession_info', {}).get('initial_release_date', 'Unknown')[:10]
    except Exception as e:
        return f"Error: {e}"

import pandas as pd
import subprocess
import os

def run_cdhit_clustering(df, sequence_column, output_prefix, cdhit_executable, identity_threshold=0.8):
    """
    Perform CD-HIT clustering on sequences in a pandas DataFrame.

    Args:
        df (pd.DataFrame): DataFrame containing sequences in a column.
        sequence_column (str): Name of the column with sequences.
        output_prefix (str): Prefix for CD-HIT output files.
        cdhit_executable (str): Path to the CD-HIT executable.
        identity_threshold (float): Sequence identity threshold for clustering (default 0.8).

    Returns:
        pd.DataFrame: DataFrame with sequences and cluster information.
    """
    # Create a FASTA file from the sequence column
    input_fasta = f"{output_prefix}_input.fasta"
    output_fasta = f"{output_prefix}_output.fasta"
    cluster_file = f"{output_prefix}_output.fasta.clstr"

    with open(input_fasta, "w") as fasta_file:
        for idx, sequence in enumerate(df[sequence_column], start=1):
            fasta_file.write(f">seq{idx}\n{sequence}\n")

    # Run CD-HIT
    cdhit_command = [
        cdhit_executable,
        "-i", input_fasta,
        "-o", output_fasta,
        "-c", str(identity_threshold),
        "-T", "0",  # Use all available threads
        "-M", "2000",  # Set memory limit to 2000 MB (greater than 1573 MB)
        "-l", "5"
    ]

    subprocess.run(cdhit_command, check=True)

    # Parse the clustering results
    clusters = {}
    current_cluster = None

    with open(cluster_file, "r") as clstr:
        for line in clstr:
            if line.startswith(">Cluster"):
                current_cluster = int(line.split()[1])
            #elif "*" in line:  # Representative sequence
            elif ">" in line:  # Sequence line
                seq_id = line.split(">")[1].split("...")[0]
                clusters[seq_id] = current_cluster

    # Map cluster information back to the DataFrame
    sequence_to_cluster = {f"seq{idx + 1}": clusters.get(f"seq{idx + 1}", None) for idx in range(len(df))}
    df["cluster"] = df.index.map(lambda idx: sequence_to_cluster.get(f"seq{idx + 1}"))

    # Clean up intermediate files
    # os.remove(input_fasta)
    # os.remove(output_fasta)
    # os.remove(cluster_file)

    return df