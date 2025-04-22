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
    #modified_to_unmodified =  
    # Manually curated mapping of noncanonical residue names to RNA.
    modified_to_unmodified = {   'A'  :'  A',   'C':'  C',   'G':'  G',   'U':'  U',\
                    '5BU':'  U', 'OMC':'  C', '5MC':'  C', 'CCC':'  C', ' DC':'  C', \
                    'CBR':'  C', 'CBV':'  C', 'CB2':'  C', '2MG':'  G', \
                    'H2U':'  U', 'PSU':'  U', '  U':'  U', '5MU':'  U', '2MU':'  U', 'OMU':'  U', \
                    'OMG':'  G', '7MG':'  G', '1MG':'  G', 'GTP':'  G', 'AMP':'  A', 'MIA':'  A', ' YG':'  G', \
                    'M2G':'  G', 'YYG':'  G', ' DG':'  G', 'G46':'  G', ' IC':'  C', ' IG':'  G',  \
                    'ZMP':'ZMP', 'YYG':'  G', '2MG':'  G', 'H2U':'  U', 'AG9':'  C', ' IU':'  U', '3TD':'  U', \
                    'A2M':'  A', '1MA':'  A', 'MA6':'  A', 'QUO':'  G', '6MZ':'  A', \
                    'C4J':'  C', '4OC':'  C', 'G7M':'  G', 'T6A':'  A', 'AET':'  A', 'I4U':'  U', 'UR3':'  U', \
                    'P7G':'  G', 'B9B':'  G', 'B8H':'  U', 'E6G':'  G', 'B8W':'  G', 'B8N':'  U', '4SU':'  U', \
                    'LV2':'  C', '4AC':'  C', 'UY4':'  A', 'I2T':'  C', '7SN':'  G', 'SUR':'  U', '7S3':'  G', \
                    'LHH':'  C', 'FHU':'  U', 'B9H':'  C', 'M1Y':'  U', 'B8Q':'  C',\
                    'M7A':'  A', 'B8K':'  G', '2PR':'  G', 'LCG':'  G', 'UFT':'  U', 'CFZ':'  C', \
                    'DA' :'  A',  'DC':'  C',  'DG':'  G',  'DT':'  U', '3AU':'  U', '9QV':'  U',' DU':'  U',
                    'CFL':'  C', 'T2T':'  T', 'N'  :'  A', 'I'  :'  G', 'GRB':'  G', 'E3C':'  C', \
                    'MMX':'  C', '1W5':'  C', '8AZ':'  G', 'B8T':'  C', 'UY1':'  U', '75B':'  U', \
                    '4DU':'  A', '5HM':'  C', '6FC':'  C', 'E7G':'  G', 'MHG':'  G', 'DU' :'  U', \
                    '56B':'  G', 'P5P':'  A', 'UMS':'  U', 'PYO':'  U', 'JMC':'  C', 'ZJS':'  A', \
                    '6IA':'  A', 'CM0':'  U', '2MA':'  A', 'RSP':'  U', 'UD5':'  U', 'MUM':'  U', \
                    'IU' :'  U', '12A':'  A', '70U':'  U', 'U8U':'  U',  'YG':'  G', 'BRU':'  U', \
                    'ATP':'  A', 'CTP':'  C', 'UTP':'  U', '5IU':'  I', 'GDP':'  G', '5IC':'  C', \
                }

    for key in modified_to_unmodified:
        modified_to_unmodified[key] = modified_to_unmodified[key].strip()
    
    # Initialize the parser and parse the CIF file
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("structure", input_cif)

    output_files = []
    ID = input_cif.split('/')[-1].split('.')[0]

    for model in structure:
        for chain in model:
            # Replace modified nucleotides with their unmodified counterparts
            to_delete=[]
            for residue in chain.get_residues():
                if residue.resname in modified_to_unmodified:
                    residue.resname = modified_to_unmodified[residue.resname]


                if residue.resname not in {"A", "U", "G", "C"}:
                    to_delete.append(residue.id)
            # print(chain.child_dict.keys())
            # exit()
            #print(to_delete)
            # Remove unwanted residues
            for res_id in to_delete:
                #print(chain.child_dict[res_id])
                chain.detach_child(res_id)
                #del chain.child_dict[res_id]

            #print(list(chain.child_dict.keys()))
            # Check if the chain contains only RNA residues
            is_rna = all(residue.resname in ["A", "C", "G", "U"] for residue in chain.get_residues())
            print([residue.resname for residue in chain.get_residues()])
            print(chain)
            print(is_rna)
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
    input_dir = "rna_structures/"
    output_dir = "split_chains"
    cif_files = glob(f"{input_dir}/4E8K*.cif")

    # Use multiprocessing to process CIF files in parallel
    cpu_count = os.cpu_count()
    print(f"Processing {len(cif_files)} CIF files using {cpu_count} processes...")
    with Pool(cpu_count) as pool:
        results = list(tqdm(pool.imap(process_cif_file, cif_files), total=len(cif_files), desc="Processing CIF files"))

    # Flatten the list of results
    output_files = [file for sublist in results for file in sublist]
    #print("Generated files:", output_files)