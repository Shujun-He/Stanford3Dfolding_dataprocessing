from glob import glob
from tqdm import tqdm
import numpy as np
import pandas as pd
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
import matplotlib.pyplot as plt
from collections import defaultdict
import pandas as pd
from utils import *

cif_files=glob("split_chains/*.cif")
# len(cif_files)

# stats=pd.read_csv("extracted_structures.csv")
# stats=stats.set_index('cif')

# structure_cutoff=0.1

# filtered_cif_files=[f for f in cif_files if f in stats.index and stats.loc[f,'structuredness']>structure_cutoff]

filtered_cif_files=cif_files
# print("after filtering by structuredness there are:", len(filtered_cif_files)) 
# print(len(filtered_cif_files))

#exit()

def get_nan_vector():
    # Create a 1x3 vector of NaNs
    nan_vector = np.full((3), np.nan)

    # Alternatively, using numpy.array
    #nan_vector_alt = np.array([np.nan, np.nan, np.nan]).reshape(1, 3)
    return nan_vector

from Bio import SeqIO

def get_sequence(file):

    # Step 2: Replace "RNAsolo_member_cif" with "solo_member_fasta"
    modified_string = file.replace(".cif", ".fasta")


    # Read the raw sequence
    with open(modified_string, "r") as file:
        # Skip the first line (header)
        file.readline()

        # Read the remaining lines and join them into a single sequence string
        sequence = ''.join(line.strip() for line in file)
        
    return sequence

rna_atom_groups = {
    "A": {
        "sugar_ring": ["C1'",],
    },
    "U": {
        "sugar_ring": ["C1'"],
    },
    "G": {
        "sugar_ring": ["C1'"],
    },
    "C": {
        "sugar_ring": ["C1'"],
    }
}


def get_xyz_sequence(file):
    pdb_info = MMCIF2Dict(file)

    xyz = [
        pdb_info['_atom_site.Cartn_x'],
        pdb_info['_atom_site.Cartn_y'],
        pdb_info['_atom_site.Cartn_z']
    ]

    xyz = np.array(xyz, dtype='float32').T

    seq_id = np.array([float(x) if x != '.' else np.nan for x in pdb_info['_atom_site.label_seq_id']], dtype='float32')
    atom_id = np.array(pdb_info['_atom_site.label_atom_id'])
    res_id = np.array(pdb_info['_atom_site.label_comp_id'])

    unique_seq_id = np.unique(seq_id[seq_id==seq_id]).astype('int')
    sequence_res = get_sequence(file)
    sequence_complete = ['N'] * int(unique_seq_id.max().item())

    if len(unique_seq_id)!=len(sequence_res): #skip some special cases like HEATOM
        return None, None, None
    
    for i, nt in zip(unique_seq_id, sequence_res):
        if i > 0:
            sequence_complete[i - 1] = nt

    grouped_xyz= []


    for i in range(len(sequence_complete)):
        res = sequence_complete[i]
        atom_groups = rna_atom_groups[res]
        atom_coords = defaultdict(list)
        
        if res in ['A','U','G','C']:
            
            for group in atom_groups:
                
                for atom_key in atom_groups[group]:

                    if atom_key!="O3'":
                        atom_index = np.where((seq_id == (i+1)) & (atom_id == atom_key))[0]
                    else:
                        atom_index = np.where((seq_id == (i)) & (atom_id == atom_key))[0] #take O3' from previous nucleotide
                    #atom_index = np.where((seq_id == (i+1)) & (atom_id == atom_key))[0]
                    if len(atom_index) > 0:
                        atom_coords[group].append(xyz[atom_index[0]])
                    else:
                        atom_coords[group].append(get_nan_vector())
        else:
            for group in atom_groups:                
                for atom_key in atom_groups[group]:
                    atom_coords[group].append(get_nan_vector())

        for group in atom_groups:
            atom_coords[group]=np.array(atom_coords[group]).astype('float32')

        grouped_xyz.append(atom_coords)

    assert len(sequence_complete) == len(grouped_xyz)

    #xyz_3bead = np.array(xyz_3bead, dtype='float32')
    #data_3_bead.append(xyz_3bead)
    # data_sequence.append(''.join(sequence_complete))
    # data_xyz.append(grouped_xyz)

    return ''.join(sequence_complete), grouped_xyz, file

# data_xyz = []
# data_sequence = []
# for file in tqdm(filtered_cif_files):
#     get_xyz_sequence(file)


from multiprocessing import Pool, cpu_count
from functools import partial
from tqdm import tqdm
import pickle

# Function to wrap `get_xyz_sequence` for multiprocessing
def process_file(file):
    try:
        return get_xyz_sequence(file)
    except Exception as e:
        print(f"Error processing {file}: {e}")
        return None, None

#if __name__ == "__main__":
    # Use all available CPU cores
num_cores = cpu_count()

#print(f"Using {num_cores} cores for multiprocessing")

# Initialize the pool of workers
with Pool(num_cores) as pool:
    # Use tqdm to track progress
    results = list(tqdm(pool.imap(process_file, filtered_cif_files), total=len(filtered_cif_files)))

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
            data_cif_files.append(f)

print(f"Processed {len(data_sequence)} sequences and their coordinates.")

# Save results as a pickle dict
with open("raw_pdb_xyz_data.pkl", "wb+") as f:
    pickle.dump({"sequence": data_sequence, "xyz": data_xyz, "data_cif_files": data_cif_files}, f)


