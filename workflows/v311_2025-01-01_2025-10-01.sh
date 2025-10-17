#!/bin/bash
# Get the relative path to the script directory from the current script even if link or run from another directory
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
SCRIPT_DIR=$(realpath "${SCRIPT_DIR}/..")
# All of those create outputs in the current directory
# Get all RNA structures from PDB in the given date range
#
START="2025-01-01"
STOP="2025-10-01"
VERSION="v311"

# Fail if any command fails
set -e

NAME_PREFIX="${VERSION}_${START}_${STOP}"
echo "Using name prefix: ${NAME_PREFIX}"

echo '================================'
echo "Fetching data from ${START} to ${STOP}"
python "${SCRIPT_DIR}/get_pdb_data.py" --start_date ${START}  --stop_date ${STOP} --search_type polymer_type --search_term RNA --output_dir rna_structures_all --cache_dir ${SCRIPT_DIR}/output/all_cache
echo '================================'

echo "Filter out structures with sequences observed before 2025-01-01 (using precomputed cutoff file)"
echo '================================'
python "${SCRIPT_DIR}/filter_by_temporal_cutoff.py" rna_structures_all rna_structures_novel --cutoff 2025-01-01 --cutoff-file ${SCRIPT_DIR}/data/temporal_cutoffs_20251001.csv
echo '================================'

echo "Filter out ribosomal structures (all text in json and mmCIF files)"
echo '================================'
python "${SCRIPT_DIR}/filter_by_keywords.py" --is_not --keyword "ribosom" rna_structures_novel rna_structures_norib 
echo '================================'

echo "Filter by composition, keep only those with at least 50% RNA by mass, used for experiment (in constrast to observed/modeled in data)"
echo '================================'
python "${SCRIPT_DIR}/filter_by_composition.py" --threshold 0.5 rna_structures_norib rna_structures_all_assemblies 
echo '================================'

echo "Filter by assembly, keep only those with single RNA entity in first assembly"
echo '================================'
python "${SCRIPT_DIR}/filter_by_assembly.py" rna_structures_all_assemblies rna_structures 
echo '================================'

echo "Split chains into separate files"
echo '================================'
python "${SCRIPT_DIR}/split_chains.py" rna_structures split_chains 
echo '================================'

### Manual step to remove problematic chains
IDS_TO_REMOVE=(
    "8XTQ_A" # Misaligned with other entry
)
for ID in "${IDS_TO_REMOVE[@]}"; do
    FILES_TO_REMOVE=$(find split_chains -type f -name "${ID}.*")
    if ls $FILES_TO_REMOVE 1> /dev/null 2>&1; then
        echo "Removing problematic files: $FILES_TO_REMOVE"
        rm $FILES_TO_REMOVE
    else
        echo "File not found, skipping removal: $FILE_TO_REMOVE"
    fi
done


echo '================================'
echo "Convert PDB files and use auth IDs (as saved by previous step using BioPython- this will match seq numbering from original files)"
echo '================================'
python "${SCRIPT_DIR}/get_xyz_data.py" split_chains raw_pdb_xyz_data.pkl --id-source auth  
echo '================================'

echo "Filter out small fragments"
echo '================================'
python "${SCRIPT_DIR}/filter_by_interC1.py" 
echo '================================'

echo "Group by sequence identity, match shorter/longer sequences if they match at least 90%"
echo '================================'
python "${SCRIPT_DIR}/groupby_sequence.py" --min_match_ratio=0.9 filtered_pdb_xyz_data.pkl grouped 
echo '================================'

echo "Generate training data keeping up to 40 structures per sequence group, using metadata from all RNA structures (to get release dates, descriptions and refernce fasta files)"
echo '================================'
python "${SCRIPT_DIR}/generate_train_multisolution_long.py" grouped test_noalign_multi_long --keep 40 --metadata-file rna_structures_all/rna_structures_metadata.json 
echo '================================'

echo "Convert to appropriate wide formats (each model/conformation in separate columns), C1' only, pad to 40 columns"
echo '================================'
python "${SCRIPT_DIR}/convert_to_multisolution_wide.py" test_noalign_multi_long_multisolution_allatom_long.csv "${NAME_PREFIX}"_C1prime.csv --select "C1'" --add-column Usage Public --drop-column group_id --fill-coordinates="-1e18" --pad 40
echo '================================'

echo "Convert to appropriate wide formats (each model/conformation in separate columns)"
echo '================================'
python "${SCRIPT_DIR}/convert_to_multisolution_wide.py" test_noalign_multi_long_multisolution_allatom_long.csv "${NAME_PREFIX}"_allatom.csv 
echo '================================'

echo "One sequence per group, keep the group names as target IDs"
echo '================================'
python "${SCRIPT_DIR}/get_single_target_per_group.py" test_noalign_multi_long_sequences.csv "${NAME_PREFIX}"_sequences.csv --sort_by temporal_cutoff target_id --rename target_id all_pdb_ids --rename group_id target_id --merge target_id description 
echo '================================'

echo "Create combined pickle file for C1' data and sequences"
echo '================================'
python "${SCRIPT_DIR}/to_pickle.py" "${NAME_PREFIX}"_C1prime.csv "${NAME_PREFIX}"_sequences.csv "${NAME_PREFIX}".pkl
echo '================================'
echo "Done, generated files with prefix ${NAME_PREFIX}_*"