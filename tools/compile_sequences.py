from collections import defaultdict
import json
import pandas as pd
from pathlib import Path
from tqdm import tqdm
import argparse

from groupby_sequence import PartialSequenceMapping


def read_fasta(fasta_file):
    """Read FASTA file and return list of sequences with headers."""
    sequences = []
    if not fasta_file or not Path(fasta_file).exists():
        return sequences

    with open(fasta_file, "r") as f:
        current_header = None
        current_sequence = []

        for line in f:
            line = line.strip()
            if line.startswith(">"):
                # Save previous sequence if exists
                if current_header is not None and current_sequence:
                    sequences.append(
                        {
                            "header": current_header,
                            "sequence": "".join(current_sequence),
                        }
                    )
                # Start new sequence
                current_header = line[1:]  # Remove '>' character
                current_sequence = []
            else:
                current_sequence.append(line)

        # Don't forget the last sequence
        if current_header is not None and current_sequence:
            sequences.append(
                {"header": current_header, "sequence": "".join(current_sequence)}
            )

    return sequences


def process_metadata_and_sequences(metadata_file, base_dir):
    """Process metadata file and compile all sequences with temporal cutoff."""

    # Load metadata
    with open(metadata_file, "r") as f:
        metadata = json.load(f)

    all_sequences = []

    print(f"Processing {len(metadata)} entries from metadata file")

    for entry in tqdm(metadata):
        pdb_id = entry.get("PDB_ID", "")
        release_date = entry.get("Release", "")
        description = entry.get("Description", "")
        title = entry.get("Title", "")
        keywords = entry.get("Keywords", "")
        fasta_file = entry.get("FASTA_File", "")

        # Convert relative path to absolute path
        if fasta_file:
            fasta_path = Path(base_dir) / fasta_file
        else:
            print(f"Warning: No FASTA file specified for {pdb_id}")
            continue

        # Read sequences from FASTA file
        sequences = read_fasta(fasta_path)

        if not sequences:
            print(f"Warning: No sequences found in {fasta_path} for {pdb_id}")
            continue

        # Process each sequence in the FASTA file
        for seq_idx, seq_data in enumerate(sequences):
            if len(set(seq_data["sequence"]) - set(["A", "C", "G", "U", "X", "N"])) > 0:
                # Skip non-RNA sequences
                continue
            # Parse header to extract chain and other info
            header_parts = seq_data["header"].split("|")
            chain_info = (
                header_parts[0] if len(header_parts) > 0 else f"{pdb_id}_{seq_idx + 1}"
            )
            chain_type = header_parts[1] if len(header_parts) > 1 else ""
            chain_description = header_parts[2] if len(header_parts) > 2 else ""
            organism = header_parts[3] if len(header_parts) > 3 else ""

            # Create target_id similar to the original format
            target_id = chain_info if "_" in chain_info else f"{pdb_id}_{seq_idx + 1}"

            sequence_entry = {
                "target_id": target_id,
                "pdb_id": pdb_id,
                "chain_info": chain_info,
                "chain_type": chain_type,
                "chain_description": chain_description,
                "organism": organism,
                "sequence": seq_data["sequence"],
                "temporal_cutoff": release_date,
                "description": description,
                "title": title,
                "keywords": keywords,
                "fasta_header": seq_data["header"],
                "sequence_length": len(seq_data["sequence"]),
            }

            all_sequences.append(sequence_entry)

    # Create DataFrame and format temporal_cutoff
    df = pd.DataFrame(all_sequences)

    # Format temporal cutoff similar to generate_train_multisolution_long.py
    if "temporal_cutoff" in df.columns:
        df["temporal_cutoff"] = pd.to_datetime(
            df["temporal_cutoff"], errors="coerce"
        ).dt.strftime("%Y-%m-%d")

    return df


def groupby_sequence(data, min_match_ratio=0.9):
    grouped_data = {}
    seq_to_id = PartialSequenceMapping(min_match_ratio=min_match_ratio)
    for i in tqdm(list(range(len(data["sequence"])))):
        seq = data["sequence"][i]
        pdb_id = data["pdb_id"][i]
        if seq not in seq_to_id:
            seq_to_id[seq] = pdb_id
            ref_id = pdb_id
            grouped_data[ref_id] = defaultdict(list)
        else:
            matched_seq, ref_id = seq_to_id.get_match_item(seq)
            # We are using partial matching, make sure the new seq is recorded, keep the original ref_id
            seq_to_id[matched_seq] = ref_id
        grouped_data[ref_id]["sequence"].append(seq)
        grouped_data[ref_id]["pdb_id"].append(pdb_id)
        grouped_data[ref_id]["temporal_cutoff"].append(data["temporal_cutoff"][i])
    return grouped_data


def main():
    parser = argparse.ArgumentParser(
        description="Compile sequences from FASTA files with metadata and temporal cutoff information"
    )
    parser.add_argument(
        "metadata_file", type=str, help="JSON file with metadata for all PDB IDs"
    )
    parser.add_argument(
        "output_file", type=str, help="Output CSV file to save compiled sequences"
    )
    parser.add_argument(
        "--base-dir",
        type=str,
        default=".",
        help="Base directory for resolving relative paths in metadata (default: current directory)",
    )

    args = parser.parse_args()

    # Process metadata and compile sequences
    df = process_metadata_and_sequences(args.metadata_file, args.base_dir)

    # Save to CSV
    df.to_csv(args.output_file, index=False)

    print(f"Compiled {len(df)} sequences from {df['pdb_id'].nunique()} PDB entries")
    print(f"Results saved to {args.output_file}")

    # Print some statistics
    print("\nSequence statistics:")
    print(f"- Total sequences: {len(df)}")
    print(f"- Unique PDB IDs: {df['pdb_id'].nunique()}")
    print(f"- Average sequence length: {df['sequence_length'].mean():.1f}")
    print(
        f"- Sequence length range: {df['sequence_length'].min()} - {df['sequence_length'].max()}"
    )

    if "temporal_cutoff" in df.columns and not df["temporal_cutoff"].isna().all():
        print(
            f"- Date range: {df['temporal_cutoff'].min()} to {df['temporal_cutoff'].max()}"
        )

    # Group by sequence with partial matching
    grouped_data = groupby_sequence(df, min_match_ratio=0.9)
    print(f"\nGrouped into {len(grouped_data)} unique sequences with partial matching")
    pdb_ids = df["pdb_id"].unique()
    # For each pdb_id, find a group in which it belongs and get the earliest temporal_cutoff
    temporal_cutoffs_per_pdb = {}
    for pdb_id in tqdm(pdb_ids):
        current_cutoff = "9999-12-31"
        for group_id, group_data in grouped_data.items():
            if pdb_id in group_data["pdb_id"]:
                temporal_cutoffs = group_data["temporal_cutoff"]
                index_min = min(
                    range(len(temporal_cutoffs)), key=temporal_cutoffs.__getitem__
                )
                cutoff = temporal_cutoffs[index_min]
                if cutoff < current_cutoff:
                    current_cutoff = cutoff
                    cutoff_pdb_id = group_data["pdb_id"][index_min]
                temporal_cutoffs_per_pdb[pdb_id] = {
                    "cutoff": current_cutoff,
                    "from_pdb_id": cutoff_pdb_id,
                }
    cutoffs = pd.DataFrame.from_dict(
        temporal_cutoffs_per_pdb,
        orient="index",
    )
    cutoffs.index.name = "pdb_id"
    cutoffs.reset_index(inplace=True)
    cutoff_file = Path(args.output_file).with_name("temporal_cutoffs.csv")
    cutoffs.to_csv(cutoff_file, index=False)
    print(f"Earliest temporal cutoffs saved to {cutoff_file}")


if __name__ == "__main__":
    main()
