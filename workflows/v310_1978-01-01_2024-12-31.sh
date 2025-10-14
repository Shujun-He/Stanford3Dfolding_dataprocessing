#!/bin/bash
# Get the relative path to the script directory from the current script even if link or run from another directory
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
SCRIPT_DIR=$(realpath "${SCRIPT_DIR}/..")
# All of those create outputs in the current directory
# Get all RNA structures from PDB in the given date range
#
START="1978-01-01"
STOP="2025-10-01"
VERSION="v310"

NAME_PREFIX="${VERSION}_${START}_${STOP}"

python "${SCRIPT_DIR}/get_pdb_data.py" --start_date ${START}  --stop_date ${STOP} --search_type polymer_type --search_term RNA --output_dir rna_structures_all &&
# Filter out ribosomal structures (all text in json and mmCIF files)
python "${SCRIPT_DIR}/filter_by_keywords.py" --is_not --keyword "ribosom" rna_structures_all rna_structures_norib &&
# Filter by composition, keep only those with at least 50% RNA by mass, used for experiment (in constrast to observed/modeled in data)
python "${SCRIPT_DIR}/filter_by_composition.py" --threshold 0.5 rna_structures_norib rna_structures_all_assemblies &&
# Filter by assembly, keep only those with single RNA entity in first assembly
python "${SCRIPT_DIR}/filter_by_assembly.py" rna_structures_all_assemblies rna_structures &&
# Split chains into separate files
python "${SCRIPT_DIR}/split_chains.py" rna_structures split_chains &&
# Convert PDB files and use auth IDs (as saved by previous step using BioPython- this will match seq numbering from original files)
python "${SCRIPT_DIR}/get_xyz_data.py" split_chains raw_pdb_xyz_data.pkl --id-source auth && 
# Filter out small fragments
python "${SCRIPT_DIR}/filter_by_interC1.py" &&
# Group by sequence identity, match shorter/longer sequences if they match at least 90%
python "${SCRIPT_DIR}/groupby_sequence.py" --min_match_ratio=0.9 filtered_pdb_xyz_data.pkl grouped &&
# Generate training data keeping up to 40 structures per sequence group, using metadata from all RNA structures (to get release dates, descriptions and refernce fasta files)
python "${SCRIPT_DIR}/generate_train_multisolution_long.py" grouped test_noalign_multi_long --keep 40 --metadata-file rna_structures_all/rna_structures_metadata.json &&
# Convert to appropriate wide formats (each model/conformation in separate columns)
python "${SCRIPT_DIR}/convert_to_multisolution_wide.py" test_noalign_multi_long_multisolution_allatom_long.csv "${NAME_PREFIX}"_C1prime.csv --select "C1'" --add-column Usage Public --drop-column group_id --fill-coordinates="-1e18" &&
python "${SCRIPT_DIR}/convert_to_multisolution_wide.py" test_noalign_multi_long_multisolution_allatom_long.csv "${NAME_PREFIX}"_allatom.csv &&
# One sequence per group, keep the group names as target IDs, 
python "${SCRIPT_DIR}/get_single_target_per_group.py" test_noalign_multi_long_sequences.csv "${NAME_PREFIX}"_sequences.csv --rename target_id all_pdb_ids --rename group_id target_id --merge target_id description &&
# Create combined pickle file for C1' data and sequences
python "${SCRIPT_DIR}/to_pickle.py" "${NAME_PREFIX}"_C1prime.csv "${NAME_PREFIX}"_sequences.csv "${NAME_PREFIX}".pkl
