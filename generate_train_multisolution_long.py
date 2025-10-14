import pickle
import json
from tqdm import tqdm
import pandas as pd
import numpy as np
from rna_reference import rna_atom_groups


def process_entry(pdb_id, data, i):
    # realease_date=data['publication_date'][i]
    # description=aligned_sequences[i]['description']
    # all_sequences=aligned_sequences[i]['all_sequences']

    sequence = data["sequence"][i]
    xyz = [{} for _ in range(len(sequence))]
    for j in range(len(sequence)):
        xyz[j] = data["xyz"][i][j]["all"]

    train_data = []

    for j in range(len(sequence)):
        resname = sequence[j]
        atom_names = (
            rna_atom_groups[resname]["all"] if resname in rna_atom_groups else []
        )
        row = {"ID": f"{pdb_id}_{j + 1}", "resname": sequence[j], "resid": j + 1}

        for k, atom_name in enumerate(atom_names):
            row[f"{atom_name}_x"] = xyz[j][k][0] if k < len(xyz[j]) else np.nan
            row[f"{atom_name}_y"] = xyz[j][k][1] if k < len(xyz[j]) else np.nan
            row[f"{atom_name}_z"] = xyz[j][k][2] if k < len(xyz[j]) else np.nan
        train_data.append(row)

        # train_data.append({"ID":f"{pdb_id}_{j+1}", 'resname':sequence[j], 'resid':j+1, 'x_1':xyz[j][0], 'y_1':xyz[j][1], 'z_1':xyz[j][2]})
    train_data = pd.DataFrame(train_data)
    return train_data


def process_all(input_file, metadata_file, keep_per_group=5):
    with open(input_file, "rb") as f:
        data = pickle.load(f)

    with open(metadata_file, "rb") as f:
        metadata = json.load(f)
        metadata = {entry["PDB_ID"]: entry for entry in metadata}

    train_solution = []
    train_sequences = []
    L = len(data["sequence"])
    if keep_per_group:
        L = min(keep_per_group, L)

    for i in range(L):
        pdb_chain_id = data["pdb_ids"][i]
        train_data = process_entry(pdb_chain_id, data, i)
        train_data["solution_id"] = i + 1
        train_data["pdb_chain_id"] = pdb_chain_id
        train_solution.append(train_data)

        # Get data from metadata file and merge
        pdb_id = pdb_chain_id.split("_")[0]
        if pdb_id in metadata:
            release_date = metadata[pdb_id].get("Release", "")
            description = metadata[pdb_id].get("Description", "")
            sequence_file = metadata[pdb_id].get("FASTA_File", "")

            all_sequences = open(sequence_file).read().strip() if sequence_file else ""
        else:
            print(f"Warning: {pdb_id} not found in metadata")
            release_date = ""
            description = ""
            all_sequences = ""

        train_sequences.append(
            {
                "target_id": pdb_chain_id,
                "sequence": data["sequence"][i],
                "temporal_cutoff": release_date,
                "description": description,
                "all_sequences": all_sequences,
            }
        )

    train_solution = pd.concat(train_solution, ignore_index=True)
    train_sequences = pd.DataFrame(train_sequences)
    train_sequences["temporal_cutoff"] = pd.to_datetime(
        train_sequences["temporal_cutoff"]
    ).dt.strftime("%Y-%m-%d")
    return train_solution, train_sequences


if __name__ == "__main__":
    import argparse
    from pathlib import Path

    # from utils import parallel_process
    # from functools import partial
    import pickle

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "input_dir",
        type=str,
        help="Directory containing pickle files with data to be converted to CSV format",
    )
    parser.add_argument(
        "output_file_prefix",
        type=str,
        help="Save results to CSV files with this prefix. Two files will be created: <prefix>_sequences.csv and <prefix>_solution.csv",
    )

    parser.add_argument(
        "--metadata-file",
        required=True,
        type=str,
        help="JSON file with metadata for all PDB IDs",
    )

    parser.add_argument(
        "--keep", type=int, default=5, help="Keep this many sequences per group"
    )
    parser.add_argument(
        "--pad", action="store_true", help="Add empty elements so groups have same size"
    )

    args = parser.parse_args()
    input_dir = Path(args.input_dir)
    output_file_prefix = args.output_file_prefix

    result_solution = pd.DataFrame()
    result_sequences = pd.DataFrame()
    group_list = list(input_dir.glob("*_grouped.pkl"))
    print(f"Processing {len(group_list)} groups from {input_dir}")
    for group_file in tqdm(group_list):
        train_solution, train_sequences = process_all(
            group_file, metadata_file=args.metadata_file, keep_per_group=args.keep
        )
        group_id = group_file.stem.replace("_grouped", "")
        train_solution["group_id"] = group_id
        train_sequences["group_id"] = group_id
        result_solution = pd.concat([result_solution, train_solution])
        result_sequences = pd.concat([result_sequences, train_sequences])

    result_solution.to_csv(
        f"{output_file_prefix}_multisolution_allatom_long.csv",
        float_format="%.3f",
        index=False,
    )
    result_sequences.to_csv(f"{output_file_prefix}_sequences.csv", index=False)
