#!/usr/bin/env python3

import pandas as pd
from itertools import product

#input_file = "WGS_ID_unanalyzed_batch1_11202024.csv"
input_file = "Manifest_410_WGS_Lung_Year_2024.csv"
output_file = "formatted_manifest_2024.csv"

# Load input CSV
df = pd.read_csv(input_file)

# Drop exact duplicate rows
df = df.drop_duplicates()

# Initialize output
rows = []
unmatched = []

# Loop through each Sherlock_PID
for pid, subdf in df.groupby("Sherlock_PID"):
    normals = subdf[subdf["Attribute"] == "Normal"]["File_ID"].tolist()
    tumors = subdf[subdf["Attribute"] == "Tumor"]["File_ID"].tolist()

    if not normals or not tumors:
        unmatched.append(pid)
        continue

    # Warn on duplication
    if len(normals) > 1:
        print(f"⚠️ {pid} has {len(normals)} Normals")
    if len(tumors) > 1:
        print(f"⚠️ {pid} has {len(tumors)} Tumors")

    # Pair all Tumor-Normal combinations
    for normal, tumor in product(normals, tumors):
        rows.append({
            "Sherlock_PID": pid,
            "filename_bam_normal": normal,
            "filename_bam_tumour": tumor
        })

# Create DataFrame from all valid pairs
final = pd.DataFrame(rows)

# Add sample IDs
final['sample_id'] = ['Sample%03d' % i for i in range(1, len(final) + 1)]

# Reorder columns
final = final[['Sherlock_PID', 'sample_id', 'filename_bam_tumour', 'filename_bam_normal']]

# Save to CSV
final.to_csv(output_file, index=False)
print(f"✅ Reformatted manifest with {len(final)} Tumor-Normal pairs saved to: {output_file}")

# Log unmatched PIDs
if unmatched:
    print("⚠️ Skipped the following Sherlock_PID(s) due to missing Tumor or Normal:")
    for pid in unmatched:
        print(f" - {pid}")
