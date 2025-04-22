import pandas as pd

# Load the CSV files
df1 = pd.read_csv('train_sequences.v0.5.1.csv')
df2 = pd.read_csv('../RibonanzaNet3D_dataprocessing/train_sequences.v0.5.1.csv')

# Print shapes
print("Shape of df1:", df1.shape)
print("Shape of df2:", df2.shape)

# Find common target_ids
common_target_ids = set(df1['target_id']).intersection(set(df2['target_id']))
print("Number of common target_ids:", len(common_target_ids))

# Find uncommon target_ids
uncommon_target_ids_df1 = sorted(set(df1['target_id']) - common_target_ids)
uncommon_target_ids_df2 = sorted(set(df2['target_id']) - common_target_ids)

# Print uncommon target_ids
print("\nUncommon target_ids in df1:")
print(uncommon_target_ids_df1)

print("\nUncommon target_ids in df2:")
print(uncommon_target_ids_df2)
