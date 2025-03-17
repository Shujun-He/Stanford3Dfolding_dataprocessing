from glob import glob
from tqdm import tqdm
import numpy as np
import pandas as pd
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
import matplotlib.pyplot as plt
from collections import defaultdict
import pandas as pd
from utils import *

cif_files=glob("split_chains/*.pdb")
len(cif_files)

stats=pd.read_csv("extracted_structures.csv")
stats['pdb']=stats['cif'].apply(lambda x: x.replace(".cif",".pdb"))
stats=stats.set_index('pdb')

structure_cutoff=0.2



filtered_cif_files=[f for f in cif_files if f in stats.index and stats.loc[f,'structuredness']>0.2]
#exit()
print("after filtering by structuredness there are:", len(filtered_cif_files)) 
#print(len(filtered_cif_files))

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
    modified_string = file.replace(".pdb", ".fasta")


    # Read the raw sequence
    with open(modified_string, "r") as file:
        # Skip the first line (header)
        file.readline()

        # Read the remaining lines and join them into a single sequence string
        sequence = ''.join(line.strip() for line in file)
        
    return sequence

rna_atom_groups = {
    "A": {
        "phosphate": ["P", "OP1", "OP2", "O5'","O3'"],
        "sugar_ring": ["C1'", "C2'", "C3'", "C4'", "O4'", "C5'"],
        "base": ["N9", "N1", "C2", "N3", "C4", "C5", "C6", "N7", "C8", "N6"]
    },
    "U": {
        "phosphate": ["P", "OP1", "OP2", "O5'","O3'"],
        "sugar_ring": ["C1'", "C2'", "C3'", "C4'", "O4'", "C5'"],
        "base": ["N1", "C2", "O2", "N3", "C4", "O4", "C5", "C6"]
    },
    "G": {
        "phosphate": ["P", "OP1", "OP2", "O5'","O3'"],
        "sugar_ring": ["C1'", "C2'", "C3'", "C4'", "O4'", "C5'"],
        "base": ["N9", "N1", "C2", "N2", "N3", "C4", "C5", "C6", "O6", "N7", "C8", ]
    },
    "C": {
        "phosphate": ["P", "OP1", "OP2", "O5'","O3'"],
        "sugar_ring": ["C1'", "C2'", "C3'", "C4'", "O4'", "C5'"],
        "base": ["N1", "C2", "O2", "N3", "C4", "N4", "C5", "C6"]
    },
    "N": {
        "phosphate": ["P", "OP1", "OP2", "O5'","O3'"],
        "sugar_ring": ["C1'", "C2'", "C3'", "C4'", "O4'", "C5'"],
        "base": ["N1", "C2", "O2", "N3", "C4", "N4", "C5", "C6"]
    }
}
from Bio.PDB import PDBParser
import numpy as np

# Initialize PDB parser
parser = PDBParser(QUIET=True)

def get_xyz_sequence(file):
    """Extracts atomic coordinates and sequences from a PDB file."""
    structure = parser.get_structure("RNA", file)

    xyz = []
    seq_id = []
    atom_id = []
    res_id = []
    
    sequence_res = get_sequence(file)  # Assuming you have a function to fetch the sequence from a FASTA file

    for model in structure:
        for chain in model:
            for residue in chain:
                res_name = residue.get_resname().strip()
                
                if res_name not in ["A", "U", "G", "C"]:  # Skip non-RNA residues
                    continue

                for atom in residue:
                    xyz.append(atom.get_coord())  # Extract coordinates
                    seq_id.append(residue.get_id()[1])  # Residue sequence ID
                    atom_id.append(atom.get_name())  # Atom name
                    res_id.append(res_name)  # Residue (nucleotide) name



    # Convert to numpy arrays
    xyz = np.array(xyz, dtype='float32')
    seq_id = np.array(seq_id, dtype='int32')
    atom_id = np.array(atom_id)
    res_id = np.array(res_id)
    #print(res_id)
    unique_seq_id = np.unique(seq_id[~np.isnan(seq_id)]).astype(int)

    # if len(unique_seq_id) == 0:
    #     return None, None, None

    #sequence_complete = ['N'] * int(unique_seq_id.max().item())
    #sequence_complete = ['N'] * len(sequence_res)
    sequence_complete = list(sequence_res)
    # if len(unique_seq_id)!=len(sequence_res): #skip some special cases like HEATOM
    #     return None, None, None
    
    # for i, nt in zip(unique_seq_id, sequence_res):
    #     if i > 0:
    #         sequence_complete[i - 1] = nt

    grouped_xyz= []
    print(len(unique_seq_id))
    print(len(sequence_complete))
    # exit()

    #for i in range(len(sequence_complete)):
    grouped_xyz = [None] * len(sequence_complete)
    for i in list(unique_seq_id-1):
        #print(i)
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
        # else:
        #     for group in atom_groups:                
        #         for atom_key in atom_groups[group]:
        #             atom_coords[group].append(get_nan_vector())

            for group in atom_groups:
                atom_coords[group]=np.array(atom_coords[group]).astype('float32')

            grouped_xyz[i]=atom_coords

    #assert len(sequence_complete) == len(grouped_xyz)

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
for sequence, grouped_xyz, f in results:
    if sequence is not None and grouped_xyz is not None:
        data_sequence.append(sequence)
        data_xyz.append(grouped_xyz)
        data_cif_files.append(f)

print(f"Processed {len(data_sequence)} sequences and their coordinates.")

# Save results as a pickle file
#data = {"sequence": data_sequence, "xyz": data_xyz}


#bespoke chain break filter based on phosphate2phosphate distance
chain_break_filter=[]
for xyz in tqdm(data_xyz):
    P_xyz=[]
    for nt_xyz in xyz:
        P_xyz.append(nt_xyz['phosphate'][0])
    P_xyz=np.array(P_xyz)
    dists=np.linalg.norm(P_xyz[1:]-P_xyz[:-1], axis=1)

    if dists.max()>10:
        chain_break_filter.append(False)
    else:   
        chain_break_filter.append(True)   
data_xyz=[data_xyz[i] for i in range(len(data_xyz)) if chain_break_filter[i]]
data_sequence=[data_sequence[i] for i in range(len(data_sequence)) if chain_break_filter[i]]
data_cif_files=[data_cif_files[i] for i in range(len(data_cif_files)) if chain_break_filter[i]]
print("after filtering by chain breaks there are:", len(data_sequence))


# Use all available CPU cores
num_cores = cpu_count()

print(f"Using {num_cores} cores for multiprocessing")
pdb_ids = [f.split("/")[-1].split(".")[0].split("_")[0] for f in data_cif_files]
publication_date = [get_pdb_publication_date(pdb_id) for pdb_id in tqdm(pdb_ids)]



# Sort by publication date and drop duplicates
sorted_indices = np.argsort(publication_date)
sorted_data_sequence = [data_sequence[i] for i in sorted_indices]
sorted_data_xyz = [data_xyz[i] for i in sorted_indices]
sorted_data_cif_files = [data_cif_files[i] for i in sorted_indices]
sorted_publication_date = [publication_date[i] for i in sorted_indices]

# Use a dictionary to drop duplicates while preserving order
unique_data = {}
for seq, xyz, cif, pub_date in zip(sorted_data_sequence, sorted_data_xyz, sorted_data_cif_files, sorted_publication_date):
    if seq not in unique_data:
        unique_data[seq] = (xyz, cif, pub_date)

filtered_data_sequence = list(unique_data.keys())
filtered_data_xyz = [unique_data[seq][0] for seq in filtered_data_sequence]
filtered_data_cif_files = [unique_data[seq][1] for seq in filtered_data_sequence]
filtered_publication_date = [unique_data[seq][2] for seq in filtered_data_sequence]

print(f"Filtered to {len(filtered_data_sequence)} unique sequences.")





data = {"sequence": filtered_data_sequence, 
    "xyz": filtered_data_xyz,
    'publication_date': filtered_publication_date,
    'filtered_cif_files': filtered_data_cif_files}




df = pd.DataFrame(data)

# Path to CD-HIT executable
cdhit_executable = "./../cdhit/cd-hit"  # Adjust the path to your CD-HIT executable
output_prefix = "cdhit_results"

clustered_df = run_cdhit_clustering(df, "sequence", output_prefix, cdhit_executable)

# Drop the "xyz" column
clustered_df_without_xyz = clustered_df.drop(columns=["xyz"], errors='ignore')

# Save to CSV
clustered_df_without_xyz.to_csv("cif_3frame_data_v3_no_xyz.csv", index=False)

# Save to Parquet
clustered_df_without_xyz.to_parquet("cif_3frame_data_v3_no_xyz.parquet")


data['cluster']=clustered_df['cluster']


with open("cif_3frame_data_v3.pkl", "wb+") as f:
    pickle.dump(data, f)