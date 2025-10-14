#!/usr/bin/env python3
"""
Script to combine XYZ coordinate data and sequence data into a single pickle file.

This script reads CSV files containing XYZ coordinates and sequence information,
processes the data to group XYZ coordinates by PDB ID, and saves the combined
data structure as a pickle file.
"""

import os
import argparse
import pickle
import pandas as pd
import numpy as np
from pathlib import Path


def process_xyz_data(xyz_df):
    """
    Process XYZ DataFrame by extracting PDB IDs and grouping coordinates.

    Args:
        xyz_df (pd.DataFrame): DataFrame containing XYZ coordinates with 'ID' column

    Returns:
        list: List of numpy arrays containing grouped XYZ coordinates
    """
    # Extract PDB ID from the ID column (everything before the last underscore)
    xyz_df["pdbid"] = xyz_df["ID"].str.rsplit("_", n=1).str[0]

    # Group XYZ coordinates by PDB ID
    grouped_coords = (
        xyz_df[["pdbid", "x_1", "y_1", "z_1"]].groupby("pdbid", sort=False).agg(list)
    )

    # Convert to numpy arrays and handle special values
    xyz_arrays = []
    for i in range(len(grouped_coords)):
        arr = np.array(grouped_coords.values[i].tolist())

        # Replace extreme values with NaN
        mask = np.isclose(arr, -1e18)
        arr[mask] = float("nan")

        # Transpose the array (from the final working version in the session)
        xyz_arrays.append(arr.T)

    return xyz_arrays


def combine_data(xyz_file, seq_file, output_file):
    """
    Combine XYZ and sequence data into a single data structure and save as pickle.

    Args:
        xyz_file (str): Path to XYZ coordinates CSV file
        seq_file (str): Path to sequences CSV file
        output_file (str): Path for output pickle file
    """
    print(f"Reading XYZ data from: {xyz_file}")
    xyz_df = pd.read_csv(xyz_file)

    print(f"Reading sequence data from: {seq_file}")
    seq_df = pd.read_csv(seq_file)

    print("Processing XYZ coordinates...")
    xyz_arrays = process_xyz_data(xyz_df)

    print("Combining data...")
    result = {
        "sequence": seq_df["sequence"].to_list(),
        "temporal_cutoff": seq_df["temporal_cutoff"].to_list(),
        "description": seq_df["description"].to_list(),
        "all_sequences": seq_df["all_sequences"].to_list(),
        "xyz": xyz_arrays,
        "target_id": seq_df["target_id"].to_list(),
    }

    print(f"Saving combined data to: {output_file}")
    with open(output_file, "wb") as fh:
        pickle.dump(result, fh)

    print(f"Successfully saved {len(xyz_arrays)} structures")


def main():
    parser = argparse.ArgumentParser(
        description="Combine XYZ coordinates and sequence data into pickle file",
    )

    parser.add_argument("xyz_file", help="Path to XYZ coordinates CSV file")

    parser.add_argument("seq_file", help="Path to sequences CSV file")

    parser.add_argument("output_file", help="Output pickle file path")

    args = parser.parse_args()

    # Validate input files exist
    if not os.path.exists(args.xyz_file):
        raise FileNotFoundError(f"XYZ file not found: {args.xyz_file}")

    if not os.path.exists(args.seq_file):
        raise FileNotFoundError(f"Sequence file not found: {args.seq_file}")

    # Create output directory if needed
    output_dir = os.path.dirname(args.output_file)
    if output_dir:
        Path(output_dir).mkdir(parents=True, exist_ok=True)

    combine_data(args.xyz_file, args.seq_file, args.output_file)


if __name__ == "__main__":
    main()
