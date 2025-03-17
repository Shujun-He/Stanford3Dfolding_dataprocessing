from Bio.PDB import MMCIFParser, MMCIFIO, PDBIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO import write
import os
from glob import glob
from multiprocessing import Pool
from tqdm import tqdm

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
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Map of modified nucleotides to their unmodified counterparts
    modified_to_unmodified = {
        "PSU": "U",  # Pseudouridine to Uracil
        # Add other mappings as needed
    }

    # Initialize the parser and parse the CIF file
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("structure", input_cif)

    output_files = []
    ID = input_cif.split('/')[-1].split('.')[0]

    for model in structure:
        for chain in model:
            # Replace modified nucleotides with their unmodified counterparts
            for residue in chain.get_residues():
                if residue.resname in modified_to_unmodified:
                    residue.resname = modified_to_unmodified[residue.resname]

            # Check if the chain contains only RNA residues
            is_rna = all(residue.resname in ["A", "C", "G", "U"] for residue in chain.get_residues())

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

                # Extract the sequence and write to a FASTA file
                sequence = "".join(residue.resname for residue in chain.get_residues())
                fasta_output_path = os.path.join(output_dir, f"{ID}_{chain_id}.fasta")
                seq_record = SeqRecord(Seq(sequence), id=f"{ID}_{chain_id}", description=f"Chain {chain_id}")
                with open(fasta_output_path, "w") as fasta_file:
                    write(seq_record, fasta_file, "fasta")
                output_files.append(fasta_output_path)

    return output_files

def process_cif_file(input_cif):
    output_dir = "split_chains"
    return split_cif_by_chains(input_cif, output_dir)

if __name__ == "__main__":
    input_dir = "rna_structures"
    output_dir = "split_chains"
    cif_files = glob(f"{input_dir}/*.cif")

    # Use multiprocessing to process CIF files in parallel
    cpu_count = os.cpu_count()
    print(f"Processing {len(cif_files)} CIF files using {cpu_count} processes...")
    with Pool(cpu_count) as pool:
        results = list(tqdm(pool.imap(process_cif_file, cif_files), total=len(cif_files), desc="Processing CIF files"))

    # Flatten the list of results
    output_files = [file for sublist in results for file in sublist]
    #print("Generated files:", output_files)