#!/usr/bin/env python3
"""
This filter removes PDBs that have other dominating components than RNA.

Uses fraction of RNA mass vs other polymer components based on entity data in mmCIF
"""

from pathlib import Path
from functools import partial


def filter(mmcif_file, keyword, is_not=False):
    """Checks if keyword is in description/title/keywords fields."""
    json_metadata_file = mmcif_file.with_name(mmcif_file.stem + ".metadata.json")
    with open(json_metadata_file) as f:
        metadata = f.read()
        metadata = metadata.lower()
    has_keyword = keyword.lower() in metadata

    if is_not:
        result = not has_keyword
    else:
        result = has_keyword

    # Make sure we want to  keep
    if result:
        # Check mmCIF file as well
        with open(mmcif_file) as f:
            mmcif_data = f.read().lower()
        has_keyword = keyword.lower() in mmcif_data
        if is_not:
            result = not has_keyword
        else:
            result = has_keyword
    return result


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
        "--keyword", required=True, type=str, help="Keyword to filter by"
    )

    parser.add_argument(
        "--is_not",
        action="store_true",
        help="If set, filter out files that contain the keyword instead of keeping them",
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
        partial(filter, keyword=args.keyword, is_not=args.is_not),
        file_list,
        desc="Filtering by keywords",
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
