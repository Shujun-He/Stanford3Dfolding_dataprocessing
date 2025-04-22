import pandas as pd
from tqdm import tqdm

df= pd.read_csv('train_sequences.v0.5.1.csv')

df=df[['target_id','sequence']]

df['pdb_id']=[s[:4] for s in df['target_id']]

df.to_csv('train_sequences.simple.csv', index=False)

#get classification from pdb
import pandas as pd
import matplotlib.pyplot as plt
import requests
from collections import Counter

# Load the CSV and extract PDB IDs
df = pd.read_csv('train_sequences.v0.5.1.csv')
df = df[['target_id', 'sequence']]
df['pdb_id'] = df['target_id'].str[:4]

# Remove duplicate PDB IDs
unique_pdb_ids = df['pdb_id'].unique()

# Function to fetch classification from RCSB PDB
def fetch_classification(pdb_id):
    url = f'https://data.rcsb.org/rest/v1/core/entry/{pdb_id}'
    try:
        response = requests.get(url)
        if response.status_code == 200:
            data = response.json()
            return data.get('struct_keywords', {}).get('pdbx_keywords', 'Unknown')
    except:
        return 'Unknown'
    return 'Unknown'

# Fetch classifications
classifications = [fetch_classification(pdb_id) for pdb_id in tqdm(unique_pdb_ids)]

# Count occurrences of each classification
classification_counts = Counter(classifications)

# Create DataFrame for plotting
classification_df = pd.DataFrame.from_dict(classification_counts, orient='index', columns=['count']).sort_values(by='count', ascending=False)

classification_df.to_csv('pdb_classification_counts.csv', index=True)

# Plot the bar chart
plt.figure(figsize=(12, 6))
classification_df.plot(kind='bar', legend=False)
plt.title('PDB Classification Counts')
plt.ylabel('Count')
plt.xlabel('Classification')
plt.xticks(rotation=45, ha='right', fontsize=4)
plt.tight_layout()
plt.savefig('pdb_classification_counts.png')



