import numpy as np
import os
from arnie.utils import convert_dotbracket_to_bp_list, convert_bp_list_to_dotbracket
from glob import glob
from Bio import Align
import pandas as pd
import matplotlib.pyplot as plt
from glob import glob
from tqdm import tqdm

def ct_to_bp_list(ct_file, header_length=1):
    ''' read a ct file into a bp_list
    assumes first line header_length is title 
    and starts reading only after that
    '''
    bps = open(ct_file).readlines()
    bps = bps[header_length:]
    bp_list = []
    for bp in bps:
        bp = bp.split()[::4][:2]
        if len(bp) != 0:
            bp = [int(nt) for nt in bp]
            if bp[1] != 0 and bp[0] < bp[1]:
                bp_list.append([bp[0] - 1, bp[1] - 1])
    return bp_list


def get_ss_from_pdb(args):

    pdb,tmp_folder=args

    #os.system(f'mkdir {tmp_folder}')
    #os.chdir(tmp_folder)
    #print(os.getcwd())
    try:
        os.system(f'../x3dna-dssr -i={pdb} > /dev/null 2>&1')
        #os.system(f'./../x3dna-dssr -i={pdb}')
        sequence=open("dssr-2ndstrs.dbn",'r').read().split('\n')[1]
        structure=open("dssr-2ndstrs.dbn",'r').read().split('\n')[2]
        
        bp=ct_to_bp_list("dssr-2ndstrs.ct")


        matrix=convert_dotbracket_to_bp_list(structure,allow_pseudoknots=True)

        # List of files to remove
        files = [
            "dssr-pairs.pdb",
            "dssr-multiplets.pdb",
            "dssr-stems.pdb",
            "dssr-stacks.pdb",
            "dssr-helices.pdb",
            "dssr-hairpins.pdb",
            "dssr-atom2bases.pdb",
            "dssr-torsions.txt",
            "dssr-splays.pdb",
            "dssr-iloops.pdb",
            "dssr-bulges.pdb",
            "dssr-Aminors.pdb",
            "dssr-2ndstrs.dbn",
            "dssr-2ndstrs.ct",
            "dssr-2ndstrs.bpseq"
        ]

        #Loop through and remove each file
        for file in files:
            if os.path.isfile(file):
                os.remove(file)
                #print(f"Removed {file}")
            else:
                pass
                #print(f"{file} does not exist")

        #os.chdir("..")

        return matrix, structure, sequence, bp
    except Exception as e:
        #os.chdir("..")
        return str(e), str(e), str(e), str(e)

def calc_structuredness(s):
    return (len(s)-s.count('.'))/len(s)


# data=pd.read_parquet("../Ribonanza/2D_finetuning/david_NAKB_dec/PDB_ATOM1_train_RNA_degformer_embeddings.parquet")
# filter_no_structure=[]
# for id,dbn in zip(data['Sequence_ID'],data['dbn']):
#     if set(list(dbn))==set('.'):
#         filter_no_structure.append(False)
#     else:
#         filter_no_structure.append(True)


# data=data.loc[filter_no_structure].reset_index(drop=True)

# data['structuredness']=data['dbn'].apply(calc_structuredness)

# plt.hist(data['structuredness'])
# plt.savefig("structuredness.png")

#exit()
SS=[]
cif_filenames=glob("split_chains/*.cif")
for cif in tqdm(cif_filenames):

    matrix, structure, sequence, bp = get_ss_from_pdb([cif,'tmp'])
    SS.append(structure)

df=pd.DataFrame()
df['dbn']=SS
df['cif']=cif_filenames
df['structuredness']=df['dbn'].apply(calc_structuredness)


df.to_csv('extracted_structures.csv',index=False)