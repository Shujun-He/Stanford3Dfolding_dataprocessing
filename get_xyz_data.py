import numpy as np
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from collections import defaultdict
from rna_reference import rna_atom_groups


def get_nan_vector():
    # Create a 1x3 vector of NaNs
    nan_vector = np.full((3), np.nan)

    # Alternatively, using numpy.array
    # nan_vector_alt = np.array([np.nan, np.nan, np.nan]).reshape(1, 3)
    return nan_vector


def get_sequence(file):
    # Step 2: Replace "RNAsolo_member_cif" with "solo_member_fasta"
    modified_string = file.with_suffix(".fasta")

    # Read the raw sequence
    with open(modified_string, "r") as file:
        # Skip the first line (header)
        file.readline()

        # Read the remaining lines and join them into a single sequence string
        sequence = "".join(line.strip() for line in file)

    return sequence


def get_xyz_sequence(file, id_source="label"):
    pdb_info = MMCIF2Dict(file)

    xyz = [
        pdb_info["_atom_site.Cartn_x"],
        pdb_info["_atom_site.Cartn_y"],
        pdb_info["_atom_site.Cartn_z"],
    ]

    xyz = np.array(xyz, dtype="float32").T

    seq_id = np.array(
        [
            float(x) if x != "." else np.nan
            for x in pdb_info[f"_atom_site.{id_source}_seq_id"]
        ],
        dtype="float32",
    )
    atom_id = np.array(pdb_info["_atom_site.label_atom_id"])
    res_id = np.array(pdb_info["_atom_site.label_comp_id"])
    unique_seq_id = np.unique(seq_id[seq_id == seq_id]).astype("int")
    sequence_res = get_sequence(file)

    # The fasta file is authoritative source for the sequence
    # use it to fill in any gaps in seq_id numbering

    # Assumes that the numbering is consistent with the fasta sequence
    # If files are processed by Biopython MMCIFParser with auth_residues=False
    # then the matching numbering should be in auth_seq_id instead of label_seq_id
    # If files come from PDB then label_seq_id should be consistent with the sequence in _entity_poly (and fasta)

    # The input files also should have all hetatm removed
    sequence_complete = sequence_res

    grouped_xyz = []

    for i in range(len(sequence_complete)):
        res = sequence_complete[i]
        atom_groups = rna_atom_groups[res]
        atom_coords = defaultdict(list)

        if res in ["A", "U", "G", "C"]:
            for group in atom_groups:
                for atom_key in atom_groups[group]:
                    atom_index = np.where((seq_id == (i + 1)) & (atom_id == atom_key))[
                        0
                    ]
                    if len(atom_index) > 0:
                        atom_coords[group].append(xyz[atom_index[0]])
                    else:
                        atom_coords[group].append(get_nan_vector())
        else:
            for group in atom_groups:
                for atom_key in atom_groups[group]:
                    atom_coords[group].append(get_nan_vector())

        for group in atom_groups:
            atom_coords[group] = np.array(atom_coords[group]).astype("float32")

        grouped_xyz.append(atom_coords)

    assert len(sequence_complete) == len(grouped_xyz)
    return "".join(sequence_complete), grouped_xyz, file

# Function to wrap `get_xyz_sequence` for multiprocessing
def process_file(file, id_source="label"):
    try:
        return get_xyz_sequence(file, id_source=id_source)
    except Exception as e:
        print(f"Error processing {file}: {e}")
        print(f"Dropping {file}; Reason: xyz_processing_error")
        return None, None


if __name__ == "__main__":
    import argparse
    from pathlib import Path
    from utils import parallel_process
    from functools import partial
    import pickle

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "input_dir", type=str, help="Path to the directory with mmCIF files to process"
    )
    parser.add_argument(
        "output_file", type=str, help="Save results to this pickle file"
    )
    parser.add_argument(
        "--id-source",
        type=str,
        default="label",
        choices=["label", "auth"],
        help="Source of sequence IDs in mmCIF file (label or auth)",
    )

    args = parser.parse_args()
    id_source = args.id_source

    cif_files = list(Path(f"{args.input_dir}").glob("*.cif"))

    process_file_partial = partial(process_file, id_source=id_source)
    results = parallel_process(process_file_partial, cif_files)

    data_xyz = []
    data_sequence = []
    data_cif_files = []
    # Separate results into sequences and coordinates
    for result in results:
        if len(result) == 3:
            sequence, grouped_xyz, f = result
            if sequence is not None and grouped_xyz is not None:
                data_sequence.append(sequence)
                data_xyz.append(grouped_xyz)
                data_cif_files.append(str(f))

    print(f"Processed {len(data_sequence)} sequences and their coordinates.")

    # Save results as a pickle dict
    with open(args.output_file, "wb+") as f:
        pickle.dump(
            {
                "sequence": data_sequence,
                "xyz": data_xyz,
                "data_cif_files": data_cif_files,
            },
            f,
        )
