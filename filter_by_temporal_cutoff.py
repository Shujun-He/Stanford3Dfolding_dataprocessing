#!/usr/bin/env python3
"""
This filter removes PDBs that have other dominating components than RNA.

Uses fraction of RNA mass vs other polymer components based on entity data in mmCIF
"""

from pathlib import Path
from functools import partial
import pandas as pd


if __name__ == "__main__":
    from utils import parallel_process
    import argparse

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
        "--cutoff-file", required=True, type=str, help="Path to the cutoff file"
    )
    parser.add_argument(
        "--cutoff", required=True, type=str, help="Cutoff date for filtering"
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

    cutoffs = pd.read_csv(args.cutoff_file).set_index("pdb_id")["cutoff"].to_dict()

    keep = [args.cutoff < cutoffs.get(f.stem.upper(), "9999-12-31") for f in file_list]
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
