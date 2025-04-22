# Workflow

this is v2 of stanford 3d folding competition data processing pipeline

key changes:
grab all pdbs that have RNA anywhere in text
relaxed filters


1. get_pdb_data.py to get pdbs from the protein data bank
2. split_chains.py to get separate chains from each pdb
3. get_xyz_data.py to get xyz data in a pickle file
4. filter_by_interC1.py filter by C1' distances, 0.2 of residues have to be close to some other residue that is over 4 bases apart
5. get_pub_dates_and_dedupe.py get publication dates and deduplicate
6. align_sequence.py align sequences to those fast files and fill in missing nts in pdbs
7. generate_train.py generate competition format train data


