#!/usr/bin/env python3
"""
This filter removes PDBs that have other dominating components than RNA.

Uses fraction of RNA mass vs other polymer components based on entity data in mmCIF
"""

import os
from pathlib import Path
import sys
from collections import defaultdict
from Bio.PDB.MMCIF2Dict import MMCIF2Dict


def get_composition(mmcif_file):
    cif_dict = MMCIF2Dict(mmcif_file)

    # Get entities in the structure
    entity_ids = cif_dict.get("_entity.id", [])
    entity_formula_weight = cif_dict.get("_entity.formula_weight", [])
    entity_pdbx_number_of_molecules = cif_dict.get(
        "_entity.pdbx_number_of_molecules", []
    )
    mass_by_id = {}
    for k, weight, copies in zip(
        entity_ids, entity_formula_weight, entity_pdbx_number_of_molecules
    ):
        if weight == "?" or copies == "?":
            raise ValueError(f"Missing weight or copies for entity {k} in {mmcif_file}")
        mass_by_id[k] = float(weight) * int(copies)

    # Get entity polymer types and number of residues
    entity_poly_entity_ids = cif_dict.get("_entity_poly.entity_id", [])
    entity_poly_types = cif_dict.get("_entity_poly.type", [])

    mass_by_type = defaultdict(float)
    for i, type in zip(entity_poly_entity_ids, entity_poly_types):
        if i in mass_by_id:
            mass_by_type[type] += mass_by_id[i]
    plain_mass_by_type = dict(mass_by_type)
    print(f"{mmcif_file}; mass by type: {plain_mass_by_type}")
    return plain_mass_by_type


def filter(mmcif_file, threshold=0.5):
    """Checks if mass of the RNA component is at least `threshold` fraction of total mass."""
    try:
        composition_by_type = get_composition(mmcif_file)
    except ValueError as e:
        print(f"Warning: {e}, rejecting {mmcif_file}")
        return False
    # Allowed keywords https://mmcif.wwpdb.org/dictionaries/mmcif_std.dic/Items/_entity_poly.type.html
    rna_mass = composition_by_type.get("polyribonucleotide", 0.0)
    total_mass = sum(composition_by_type.values())
    fraction_rna = rna_mass / total_mass if total_mass > 0 else 0.0
    if fraction_rna >= threshold:
        return True
    else:
        return False


if __name__ == "__main__":
    from utils import parallel_process
    import argparse
    from functools import partial

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "input_dir", type=str, help="Path to the directory with mmCIF files to check"
    )
    parser.add_argument(
        "output_dir",
        type=str,
        help="Path to the directory where filtered files will be saved (symlinked)",
    )
    parser.add_argument(
        "--threshold", type=float, default=0.5, help="Minimum fraction of RNA mass"
    )
    parser.add_argument(
        "--copy", action="store_true", help="Copy files instead of symlinking"
    )
    parser.add_argument(
        "--extra_extensions",
        nargs="+",
        default=[".fasta", ".pdb", ".metadata.json"],
        help="Extra file extensions to look for when copying files",
    )

    args = parser.parse_args()
    input_dir = Path(args.input_dir)
    file_list = input_dir.glob(
        "*.cif"
    )  # List all mmCIF files in the directory and subdirectories
    file_list = list(file_list)
    keep = parallel_process(
        partial(filter, threshold=args.threshold),
        file_list,
        desc="Filtering by composition",
    )
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"Keeping {sum(keep)} out of {len(keep)} files")

    if args.copy:
        import shutil

        symlink_or_copy = shutil.copy
    else:

        def symlink_or_copy(src, dst):
            relative_source = src.resolve().relative_to(
                dst.parent.resolve(), walk_up=True
            )
            if dst.exists():
                dst.unlink()
            dst.symlink_to(relative_source)

    for f, k in zip(file_list, keep):
        if k:
            for extension in [".cif"] + args.extra_extensions:
                in_path = f.with_suffix(extension)
                out_path = output_dir / in_path.name
                if in_path.exists() and in_path.is_file():
                    print(f"Keeping {in_path} -> {out_path}")
                    symlink_or_copy(in_path, out_path)
