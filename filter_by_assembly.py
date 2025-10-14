#!/usr/bin/env python3
"""
This filter removes PDBs that have complex biological, multichain assembly.

Uses the _pdbx_struct_assembly category in mmCIF files.
"""

from pathlib import Path
from functools import partial
from Bio.PDB.MMCIF2Dict import MMCIF2Dict


def filter(
    mmcif_file,
):
    """Checks if keyword is in description/title/keywords fields."""
    # Check mmCIF file as well
    mmcif_dict = MMCIF2Dict(str(mmcif_file))
    assemblies = mmcif_dict.get("_pdbx_struct_assembly_gen.assembly_id", [])
    oper_expressions = mmcif_dict.get("_pdbx_struct_assembly_gen.oper_expression", [])
    assemblies_chain_ids = mmcif_dict.get("_pdbx_struct_assembly_gen.asym_id_list", [])
    #####
    # Remove if we have more than one assembly per file, or if assemblies are not defined
    if not assemblies:
        reason = "not_defined"
        return {"file": mmcif_file, "keep": False, "reason": reason}

    # Not doing it now
    # elif len(assemblies) > 1:
    #     reason = "multiple_assemblies"
    #     return {"file": mmcif_file, "keep": False, "reason": reason}
    # else:
    #     # If we have exactly one assembly, we can proceed
    #     pass

    #####
    # Reject if the assembly is created by applying symmetry operators
    if any(op != "1" for op in oper_expressions):
        reason = "symmetry_operators"
        return {"file": mmcif_file, "keep": False, "reason": reason}

    #####
    # Check if the chains are RNA chains, if there are more then one RNA chain in the assembly then reject
    chains = [set(c.split(",")) for c in assemblies_chain_ids]
    rna_chains = []

    entity_ids = mmcif_dict.get("_entity_poly.entity_id", [])
    strand_ids = mmcif_dict.get("_entity_poly.pdbx_strand_id", [])
    types = mmcif_dict.get("_entity_poly.type", [])

    rna_entity_set = set()
    # Processing first assembly only, if multiple it gets rejected earlier
    for chain_id in chains[0]:
        # Check if the chain is RNA by looking at entity_poly
        for entity_id, strand_id, type in zip(entity_ids, strand_ids, types):
            if type == "polyribonucleotide" and chain_id in strand_id.split(","):
                rna_entity_set.add(entity_id)
    if len(rna_entity_set) > 1:
        reason = "multiple_rna_entity_assembly"
        return {"file": mmcif_file, "keep": False, "reason": reason}

    return {"file": mmcif_file, "keep": True, "reason": "ok"}


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
    keep = parallel_process(filter, file_list, desc="Filtering by assembly")
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    for r in keep:
        print(f"{r['file']}: keep={r['keep']} reason={r['reason']}")

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
        if k["keep"]:
            for extension in [".cif"] + args.extra_extensions:
                in_path = f.with_suffix(extension)
                out_path = output_dir / in_path.name
                if in_path.exists() and in_path.is_file():
                    print(f"Keeping {in_path} -> {out_path}")
                    symlink_or_copy(in_path, out_path)
