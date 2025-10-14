#!/bin/bash
#SBATCH --cpus-per-task=64
#SBATCH --job-name=split-chains
#SBATCH --partition cpuq

SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"

python ${SCRIPT_DIR}/get_pdb_data.py --start_date now-1w --stop_date now --search_type polymer_type --search_term RNA --output_dir rna_structures_all &&
python ${SCRIPT_DIR}/filter_by_composition.py rna_structures_all rna_structures &&
python ${SCRIPT_DIR}/split_chains.py rna_structures split_chains &&
python ${SCRIPT_DIR}/get_xyz_data.py split_chains raw_pdb_xyz_data.pkl --id-source auth &&
python ${SCRIPT_DIR}/filter_by_interC1.py &&
python ${SCRIPT_DIR}/get_pub_dates_and_dedupe.py &&
python ${SCRIPT_DIR}/align_sequence.py &&
python ${SCRIPT_DIR}/generate_train.py