import pickle
import numpy as np
from tqdm import tqdm
import multiprocessing as mp


# Save results as a pickle dict
with open("raw_pdb_xyz_data.pkl", "rb") as f:
    data=pickle.load(f)

#save cif files to csv
import pandas as pd
df=pd.DataFrame(data['data_cif_files'])
df.to_csv("raw_pdb_cif_files.csv", index=False)


def filter_condition(i):
    xyz=data['xyz'][i]

    xyz=[nt['sugar_ring'][0] for nt in xyz]
    xyz=np.stack(xyz)
    #exit()
    #first filter if percent nan > 0.5

    if np.isnan(xyz).mean()>0.5:
        return False

    #second filter if 0.2 of nts have min inter C1 distance < 12
    distance_matrix=np.square(xyz[None,:,:]-xyz[:,None,:]).sum(-1)**0.5
    distance_matrix[np.isnan(distance_matrix)]=999999
    
    #mask diagonal and off diagonal |i-j|<4
    for i in range(len(distance_matrix)):
        distance_matrix[i,i]=999999
        for j in range(i+1, len(distance_matrix)):
            if abs(i-j)<4:
                distance_matrix[i,j]=999999
                distance_matrix[j,i]=999999


    min_distances=distance_matrix.min(0)

    if (min_distances<12).mean()>0.2:
        return True
    else:
        return False


with mp.Pool(mp.cpu_count()) as pool:
    results = list(tqdm(pool.imap(filter_condition, range(len(data['xyz']))), total=len(data['xyz'])))

# Indices that passed the filter (True = keep)
passed_indices = [i for i, keep in enumerate(results) if keep]

# Optionally save or process these indices
print(f"{len(passed_indices)} entries passed the filter.")
    

#filter data
data_xyz=[data['xyz'][i] for i in passed_indices]
data_sequence=[data['sequence'][i] for i in passed_indices]
data_cif_files=[data['data_cif_files'][i] for i in passed_indices]
print("after filtering inter C1' distances there are:", len(data_sequence))

# Save results as a pickle dict
with open("filtered_pdb_xyz_data.pkl", "wb+") as f:
    pickle.dump({"sequence": data_sequence, "xyz": data_xyz, "data_cif_files": data_cif_files}, f)

# Save cif files to csv
df=pd.DataFrame(data_cif_files)
df.to_csv("filtered_pdb_cif_files.csv", index=False)