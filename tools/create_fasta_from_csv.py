#!/usr/bin/env python3
import pandas as pd
import sys
from pathlib import Path

# Configuration
CSV_PATH = Path(sys.argv[1])
OUTPUT_DIR = Path(sys.argv[2])

# Create output directory
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Read CSV
print(f"Reading CSV from {CSV_PATH}")
df = pd.read_csv(CSV_PATH)

# Verify columns
if "target_id" not in df.columns or "sequence" not in df.columns:
    print("ERROR: CSV must contain 'target_id' and 'sequence' columns")
    sys.exit(1)

print(f"Found {len(df)} sequences")

# Generate FASTA files
for idx, row in df.iterrows():
    target_id = row["target_id"]
    sequence = row["sequence"]

    # Use target_id as filename
    fasta_file = OUTPUT_DIR / f"{target_id}.fasta"

    with open(fasta_file, "w") as f:
        f.write(f">{target_id}\n")
        f.write(f"{sequence}\n")

    print(f"Created {fasta_file.name}")

print(f"\nGenerated {len(df)} FASTA files in {OUTPUT_DIR}")
