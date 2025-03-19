# Workflow

1. get_pdb_data.py to get pdbs from the protein data bank
2. split_chains.py to get separate chains from each pdb
3. get_xyz_data.py to get xyz data in a pickle file
4. align_sequence.py to align sequences to those fast files and fill in missing nts in pdbs
5. generate_train.py to generate kaggle format train data

I've uploaded the pickle file: https://www.kaggle.com/datasets/shujun717/stanford3d-dataprocessing-pickle, so you may start from step 4 to reproduce Stanford3D competition data