import copy
from Bio.PDB import MMCIFParser, MMCIFIO, PDBIO
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO import write
import os
from glob import glob
import numpy as np
import re
from rna_reference import modified_to_unmodified_nakb


def get_sequence_from_cif(cif_dict):
    chain_to_sequence = {}
    # Get sequence from entity_poly
    entity_ids = cif_dict.get("_entity_poly.entity_id", [])
    # Using canonical one-letter code that maps modified residues to parent according to wwPDB components.cif
    seqs = cif_dict.get("_entity_poly.pdbx_seq_one_letter_code_can", [])
    seqs_full = cif_dict.get("_entity_poly.pdbx_seq_one_letter_code", [])

    strand_ids = cif_dict.get("_entity_poly.pdbx_strand_id", [])
    types = cif_dict.get("_entity_poly.type", [])

    for entity_id, seq, seq_full, strand_ids, type in zip(
        entity_ids, seqs, seqs_full, strand_ids, types
    ):
        print(f"Entity {entity_id} type {type} strand_ids {strand_ids} seq {seq}")
        if type != "polyribonucleotide":
            continue
        strand_ids = strand_ids.strip().split(",")
        seq = seq.replace("\n", "").replace(" ", "").strip()

        # Check if all residues in the sequence are valid RNA residues, if other than AUCG then use seqs_full to get the ligand id and map using nakb_modified_nt
        if not set(seq) == {"A", "U", "G", "C"}:
            print(f"Entity {entity_id} has unresolved modified residues.")
            # The seq_full has modified residue names in parantheses, e.g. A(5MC)UCG
            # Extract the residue names and map them using modified_to_unmodified_nakb

            residues = re.findall(r"([A-Z]|\(.*?\))", seq_full)
            residues = [
                res[1:-1] if res.startswith("(") and res.endswith(")") else res
                for res in residues
            ]
            # Map the residues to unmodified using nakb mapping
            # If not found in the mapping use 'A'
            seq_mapped = []
            for res in residues:
                if res in {"A", "U", "G", "C"}:
                    seq_mapped.append(res)
                elif res in modified_to_unmodified_nakb:
                    print(
                        f"Mapping modified residue {res} to {modified_to_unmodified_nakb[res]} using NAKB data"
                    )
                    seq_mapped.append(modified_to_unmodified_nakb[res])
                else:
                    print(f"Warning: residue {res} not in NAKB mapping")
                    seq_mapped.append(res)
            seq = "".join(seq_mapped)

        for strand_id in strand_ids:
            if strand_id in chain_to_sequence:
                print(
                    f"Warning: multiple sequences for chain {strand_id} in {input_cif}"
                )
            chain_to_sequence[strand_id] = seq
    return chain_to_sequence


def split_cif_by_chains(input_cif, output_dir):
    """
    Splits a CIF file containing multiple chains into separate CIF files, saving only RNA chains,
    and writes the sequences of each chain to a FASTA file.

    Parameters:
        input_cif (str): Path to the input CIF file.
        output_dir (str): Directory to save the individual chain CIF files and FASTA files.

    Returns:
        list: Paths to the generated CIF and FASTA files.
    """

    # Initialize the parser and parse the CIF file
    print("Processing", input_cif)
    parser = MMCIFParser(QUIET=True, auth_residues=False)
    cif_dict = MMCIF2Dict(input_cif)
    chain_to_sequence = get_sequence_from_cif(cif_dict)
    structure = parser.get_structure("structure", input_cif)

    output_files = []
    ID = input_cif.name.split(".")[0]  # to take care of .cif.gz too

    for model in structure:
        for chain in model:
            if chain.id not in chain_to_sequence:
                # Skip non-RNA chains
                print("Available chains:", list(chain_to_sequence.keys()))
                print(f"Skipping non-RNA chain {chain.id} in {input_cif}")
                continue
            else:
                current_sequence = chain_to_sequence[chain.id]
                output_sequence = list(copy.copy(current_sequence))

            # Process residues in the chain, and assign modified residues to unmodified
            # The due to auth_residues=False in MMCIFParser, the residues should have numbering
            # consistent with the sequence in _entity_poly
            # If not raises error. If the residues is modified with unknown parent (X) then
            # use custom modified_to_unmodified to map to unmodified

            to_delete=[]
            for residue in chain.get_residues():
                # Get the single code based on the resid
                single_code = current_sequence[residue.id[1] - 1]
                # Should be A,U,G,C at this point
                if single_code not in {"A", "U", "G", "C"}:
                    print(
                        f"Warning: residue {residue.resname} at position {residue.id[1]} in chain {chain.id} not mapped to standard base"
                    )
                # TODO: verify atom mappings

                residue.resname = single_code

                # print(residue.type)

                # if residue.id[0] != " ":
                #     print(residue.resname)
                # exclude heteroatoms and non canonical residues
                # if residue.resname not in {"A", "U", "G", "C"} or residue.id[0] != " ":
                #     to_delete.append(residue.id)
            #exit()
            # print(chain.child_dict.keys())
            # exit()
            #print(to_delete)
            # Remove unwanted residues
            for res_id in to_delete:
                #print(chain.child_dict[res_id])
                chain.detach_child(res_id)
                #del chain.child_dict[res_id]

            #print(list(chain.child_dict.keys()))
            # FIXME: configurable minimum length 
            # continue if the chain has at least 10 residues left
            if len(chain.child_dict) < 10:
                print(f"Dropping {ID}_{chain_id}; Reason: rna_chain_too_short")
                continue

            # Check if 0.5 or more of the residues an the chain are RNA
            # FIXME: configurable threshold
            # TODO:idealy that should depend on the interactions, ie. it is fine to keep 
            #      RNA from protein complex if the RNA is interacting in limited fashion with the protein
            #      alternatively assign score for further filtering
            # FIXME: verify, these should be only A,U,G,C now after mapping and deleting others
            is_rna = np.array([residue.resname in ["A", "C", "G", "U"] for residue in chain.get_residues()]).mean() >= 0.5
            # print([residue.resname for residue in chain.get_residues()])  
            # print(chain)
            # print(is_rna)
            if is_rna:
                chain_id = chain.id

                # Create a new structure object with only this chain
                chain_structure = structure.__class__("{}_chain_{}".format(structure.id, chain_id))
                model_copy = model.__class__(model.id)
                model_copy.add(chain.copy())
                chain_structure.add(model_copy)

                # Write the chain to a new CIF file
                cif_output_path = os.path.join(output_dir, f"{ID}_{chain_id}.cif")
                io = MMCIFIO()
                io.set_structure(chain_structure)
                io.save(cif_output_path)
                output_files.append(cif_output_path)

                # Write the chain to a new PDB file
                try:
                    pdb_output_path = os.path.join(output_dir, f"{ID}_{chain_id}.pdb")
                    pdb_io = PDBIO()
                    pdb_io.set_structure(chain_structure)
                    pdb_io.save(pdb_output_path)
                    output_files.append(pdb_output_path)
                except:
                    pass

                # Use the modified sequence to a FASTA file
                sequence = "".join(output_sequence)
                fasta_output_path = os.path.join(output_dir, f"{ID}_{chain_id}.fasta")
                seq_record = SeqRecord(Seq(sequence), id=f"{ID}_{chain_id}", description=f"Chain {chain_id}")
                with open(fasta_output_path, "w") as fasta_file:
                    write(seq_record, fasta_file, "fasta")
                output_files.append(fasta_output_path)
                # TODO: verify that this sequence matches the one in downloaded FASTA file
            else:
                print(f"Dropping {ID}_{chain_id}; Reason: rna_chain_fraction_below_0.5")

    return output_files

if __name__ == "__main__":
    import argparse
    from functools import partial
    from pathlib import Path
    from utils import parallel_process

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "input_dir", type=str, help="Path to the directory with mmCIF files to check"
    )
    parser.add_argument(
        "output_dir",
        type=str,
        help="Path to the directory where filtered files will be saved (symlinked)",
    )
    args = parser.parse_args()
    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    cif_files = list(input_dir.glob("*.cif"))
    # Use multiprocessing to process CIF files in parallel
    parallel_process(
        partial(split_cif_by_chains, output_dir=output_dir),
        cif_files,
        desc="Splitting CIF files by chains",
    )