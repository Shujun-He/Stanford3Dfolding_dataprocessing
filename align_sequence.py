import pickle
import os
import requests
from Bio import SeqIO
from Bio import pairwise2
from tqdm import tqdm

# Load the pickle file
with open('pdb_xyz_data.pkl', 'rb') as f:
    data = pickle.load(f)

# Ensure the output directory exists
os.makedirs("alignments", exist_ok=True)

# Base URL for RCSB PDB Data API
pdb_api_url = "https://data.rcsb.org/rest/v1/core/entry/"

cnt=0
# Iterate through each sequence and corresponding CIF file
aligned_sequences = []
for sequence, file in tqdm(zip(data['sequence'], data['filtered_cif_files']), total=len(data['sequence'])):
    file_name = os.path.basename(file)
    pdb_id = file_name.split('_')[0]
    chain_id = file_name.split('_')[1].strip(".pdb")

    fasta_file_path = f"../RibonanzaNet3D_dataprocessing/rna_structures/{pdb_id}.fasta"
    alignment_file_path = f"alignments/{pdb_id}_{chain_id}_best_alignment.txt"

    # Fetch PDB entry details
    try:
        response = requests.get(f"{pdb_api_url}{pdb_id}")
        response.raise_for_status()
        pdb_data = response.json()
        description = pdb_data.get('struct', {}).get('title', 'No description available.')
    except requests.exceptions.RequestException as e:
        print(f"Error fetching data for PDB ID {pdb_id}: {e}")
        description = 'No description available.'

    best_score = -float("inf")
    best_fasta_sequence = None
    best_alignment = None
    index_map = []

    # Check if the fasta file exists
    if not os.path.exists(fasta_file_path):
        print(f"FASTA file {fasta_file_path} not found. Skipping...")
        #continue
        fasta_string = ''


    else:
        fasta_string = open(fasta_file_path).read()
        # Read sequences from FASTA file
        for record in SeqIO.parse(fasta_file_path, "fasta"):
            fasta_sequence = str(record.seq)

            #if len(fasta_sequence) < 1000:
            if len(fasta_sequence) < 1000 and len(sequence) < 1000:
                # Align with Smith-Waterman (local alignment)
                alignments = pairwise2.align.localms(sequence, fasta_sequence, 
                                                    match=1, mismatch=-float("inf"), 
                                                    open=-1, extend=-0.5)

                # Get the best alignment
                if alignments and alignments[0].score > best_score:
                    best_score = alignments[0].score
                    best_fasta_sequence = fasta_sequence
                    best_alignment = alignments[0]

        # Save the best alignment to a text file
        
        if best_alignment:
            aligned_seq1, aligned_seq2, score, start, end = best_alignment

            # Compute index mapping
            
            seq_idx = 0  # Index in `sequence`
            fasta_idx = 0  # Index in `best_fasta_sequence`
            
            for s1, s2 in zip(aligned_seq1, aligned_seq2):
                if s1 != '-':  # s1 is from `sequence`
                    if s2 != '-':  # s2 is from `best_fasta_sequence`
                        index_map.append((seq_idx, fasta_idx))
                    seq_idx += 1
                if s2 != '-':  # s2 is from `best_fasta_sequence`
                    fasta_idx += 1

            with open(alignment_file_path, "w") as f:
                f.write(f"Best alignment for PDB: {pdb_id}, Chain: {chain_id}\n")
                f.write(f"Description: {description}\n\n")
                f.write(pairwise2.format_alignment(*best_alignment))
                f.write("\n\nIndex Mapping (original sequence to aligned sequence indices):\n")
                f.write(str(index_map))

    
    aligned_sequences.append({
        "pdb_id": pdb_id,
        "chain_id": chain_id,
        "description": description,
        "best_alignment": best_alignment,
        "index_map": index_map,
        "all_sequences": fasta_string
    })
    #exit()
    cnt+=1
    # if cnt>20:
    #     exit()
#save aligned sequences
with open('aligned_sequences.pkl', 'wb') as f:
    pickle.dump(aligned_sequences, f)