import pickle
import os
import requests
from Bio import SeqIO
from Bio import pairwise2
from tqdm import tqdm
import pandas as pd
import numpy as np

# Load the pickle file
with open('deduped_pdb_xyz_data.pkl', 'rb') as f:
    data = pickle.load(f)


# for f in data['filtered_cif_files']:
#     if '7M57' in f:
#         print(f)
#         exit()
# exit()

#load aligned sequences
with open('aligned_sequences.pkl', 'rb') as f:
    aligned_sequences = pickle.load(f)


wrong_alignments = 0
train_sequences = []
train_solution = []
for i in tqdm(range(len(aligned_sequences))):
    alignment=aligned_sequences[i]['best_alignment']


    #pdb_id=aligned_sequences[i]['pdb_id']
    pdb_id=data['filtered_cif_files'][i].split('/')[1].split('.')[0]#.strip(".cif")
    #exit()
    realease_date=data['publication_date'][i]
    description=aligned_sequences[i]['description']
    all_sequences=aligned_sequences[i]['all_sequences']
    sequence=data['sequence'][i]

    if alignment is None:
        wrong_alignments+=1
        alignment_check=False
    
    elif alignment is not None:
        alignment_check = True
        # if len(alignment[1]) < len(alignment[0]):
        if alignment[0]!=alignment[1]:
            seqA = alignment[0]
            seqB = alignment[1]
            #check non-gapped positions are the same
            for j in range(len(seqA)):
                if seqA[j] != '-' and seqB[j] != '-':
                    if seqA[j] != seqB[j]:
                        alignment_check=False
                        wrong_alignments+=1
                        break
        
    if alignment_check:
        #pdb_sequence=alignment[0]
        #ref_sequence=alignment[1]

        non_gapped_indices=np.where(np.array(list(alignment[0]))!='-')[0]
        assert len(non_gapped_indices)==len(data['xyz'][i])

        # len(sequence)x3 nan default array
        xyz=np.full((len(alignment[0]),3), np.nan)
        for j,k in enumerate(non_gapped_indices):
            xyz[k]=data['xyz'][i][j]['sugar_ring'][0]
        sequence=alignment[1]
    else:
        xyz=np.full((len(sequence),3), np.nan)
        for j in range(len(sequence)):
            xyz[j]=data['xyz'][i][j]['sugar_ring'][0]

    train_data=[]
    for j in range(len(sequence)):
        train_data.append({"ID":f"{pdb_id}_{j+1}", 'resname':sequence[j], 'resid':j+1, 'x_1':xyz[j][0], 'y_1':xyz[j][1], 'z_1':xyz[j][2]})
    train_data=pd.DataFrame(train_data)
    train_solution.append(train_data)
    train_sequences.append({'target_id':pdb_id, 'sequence':sequence, 'temporal_cutoff':realease_date, 'description':description, 'all_sequences':all_sequences})

    for nt1, nt2 in zip(train_data['resname'],sequence):
        assert nt1==nt2

    #write solution df


print("Number of wrong alignments: ", wrong_alignments)
print("Used PDB sequence for all sequences with wrong alignments")

# Save the aligned sequences
train_sequences=pd.DataFrame(train_sequences)
train_solution=pd.concat(train_solution)

#write solution df
train_sequences.to_csv('train_sequences.v0.5.1.csv', index=False)
train_solution.to_csv('train_solution.v0.5.1.csv', index=False)