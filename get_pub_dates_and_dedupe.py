import pickle
import numpy as np
from tqdm import tqdm
from utils import get_pdb_publication_date
import pandas as pd

with open("filtered_pdb_xyz_data.pkl", "rb") as f:
    data=pickle.load(f)

df=pd.DataFrame(data)

# Group by 'sequence' and count how many times each appears
group_counts = df.groupby('sequence').size().reset_index(name='count')

# Sort by count descending
sorted_groups = group_counts.sort_values(by='count', ascending=False)

# Show top groups
print(sorted_groups)

#exit()

print("before filtering by unique sequences there are:", len(data['sequence']))

seen_seqs=set()
unique_indices=[]
for i in range(len(data['sequence'])):
    seq=data['sequence'][i]
    if seq not in seen_seqs:
        unique_indices.append(i)
        seen_seqs.add(seq)
#filter data
data_sequence=[data['sequence'][i] for i in unique_indices]
data_xyz=[data['xyz'][i] for i in unique_indices]
data_cif_files=[data['data_cif_files'][i] for i in unique_indices]
print("after filtering by unique sequences there are:", len(data_sequence))

#unpack
# data_sequence=data['sequence']
# data_xyz=data['xyz']
# data_cif_files=data['data_cif_files']

pdb_ids = [f.split("/")[-1].split(".")[0].split("_")[0] for f in data_cif_files]
publication_date = [get_pdb_publication_date(pdb_id) for pdb_id in tqdm(pdb_ids)]

data = {"sequence": data_sequence,
    "xyz": data_xyz,
    'publication_date': publication_date,
    'filtered_cif_files': data_cif_files}

#save
with open("deduped_pdb_xyz_data.pkl", "wb+") as f:
    pickle.dump(data, f)

#save cif files to csv
df=pd.DataFrame(data['filtered_cif_files'])
df.to_csv("deduped_pdb_xyz_data.csv", index=False)