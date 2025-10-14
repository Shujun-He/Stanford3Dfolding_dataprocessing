import gzip
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from rna_reference import modified_to_unmodified


def clean_res_name(res_name):
    striped_res_name = res_name.strip()
    if striped_res_name in ["A", "C", "G", "U"]:
        return striped_res_name
    elif striped_res_name in modified_to_unmodified:
        return modified_to_unmodified[striped_res_name].strip()
    else:  # can be modified residue with 3-letter name.
        return "X"


def extract_rna_sequence(cif_path, chain_id):
    """
    Extract RNA sequence and residue numbers from a CIF file for a specific chain.

    Modified from https://github.com/DasLab/pdb_map/blob/main/create_pdb_sequences_csv.py
    """
    if cif_path.endswith(".gz"):
        with gzip.open(cif_path, "rt") as cif_file:
            mmcif_dict = MMCIF2Dict(cif_file)
    else:
        with open(cif_path, "rt") as cif_file:
            mmcif_dict = MMCIF2Dict(cif_file)

    pdb_sequence = None
    pdb_chain_id = None
    chain_seq_nums = None
    # Extract _pdbx_poly_seq_scheme information
    strand_id = mmcif_dict.get("_pdbx_poly_seq_scheme.pdb_strand_id", [])
    mon_id = mmcif_dict.get("_pdbx_poly_seq_scheme.mon_id", [])
    pdb_mon_id = mmcif_dict.get("_pdbx_poly_seq_scheme.pdb_mon_id", [])
    pdb_seq_num = mmcif_dict.get("_pdbx_poly_seq_scheme.pdb_seq_num", [])
    auth_seq_num = mmcif_dict.get("_pdbx_poly_seq_scheme.auth_seq_num", [])
    pdb_ins_code = mmcif_dict.get("_pdbx_poly_seq_scheme.pdb_ins_code", [])
    chain_ids = list(set(strand_id))
    seq_chains = []

    full_sequence = ""
    pdb_chain_sequence = ""
    pdb_chain_seq_nums = []
    pdb_chain_ins_codes = []

    for strand, mon, pdb_mon, pdb_num, auth_num, ins_code in zip(
        strand_id, mon_id, pdb_mon_id, pdb_seq_num, auth_seq_num, pdb_ins_code
    ):
        if strand == chain_id:
            full_sequence += clean_res_name(mon)
            pdb_chain_sequence += clean_res_name(pdb_mon)
            # note use of auth_seq_num instead of pdb_seq_num since that is what Biopython uses for Residue.id
            pdb_chain_seq_nums.append(auth_num)
            pdb_chain_ins_codes.append(ins_code)

    return full_sequence, pdb_chain_sequence, pdb_chain_seq_nums, pdb_chain_ins_codes


if __name__ == "__main__":
    import sys

    cif_path = sys.argv[1]
    chain_id = sys.argv[2]
    full_sequence, pdb_chain_sequence, pdb_chain_seq_nums, pdb_chain_ins_codes = (
        extract_rna_sequence(cif_path, chain_id)
    )
    print(f"Full sequence: {full_sequence}")
    print(f"PDB chain sequence: {pdb_chain_sequence}")
    print(f"PDB chain seq nums: {pdb_chain_seq_nums}")
