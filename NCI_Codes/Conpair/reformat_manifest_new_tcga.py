#!/usr/bin/env python3

import pandas as pd
from itertools import product

input_file = "New_TCGA_Manifest_Permissions.csv"
output_file = "formatted_tcga_manifest_2024.csv"

# Load the input CSV
df = pd.read_csv(input_file)

# Drop exact duplicate rows
df = df.drop_duplicates()

# Initialize output
rows = []
unmatched = []

# Loop through each case_id
for case_id, subdf in df.groupby("cases.case_id"):
    normals = subdf[subdf["samples.tissue_type"] == "Normal"]["CRAM_File_Path"].tolist()
    tumors = subdf[subdf["samples.tissue_type"] == "Tumor"]["CRAM_File_Path"].tolist()

    if not normals or not tumors:
        unmatched.append(case_id)
        continue

    # Warn if multiple tumor/normal files found
    if len(normals) > 1:
        print(f"⚠️ {case_id} has {len(normals)} Normal samples")
    if len(tumors) > 1:
        print(f"⚠️ {case_id} has {len(tumors)} Tumor samples")

    # Pair all combinations of tumor-normal for the case
    for normal, tumor in product(normals, tumors):
        rows.append({
            "case_id": case_id,
            "filename_bam_normal": normal,
            "filename_bam_tumour": tumor
        })

# Create DataFrame
final = pd.DataFrame(rows)

# Add unique sample IDs
final["sample_id"] = ["Sample%03d" % i for i in range(1, len(final) + 1)]

# Reorder columns
final = final[["case_id", "sample_id", "filename_bam_tumour", "filename_bam_normal"]]

# Save to CSV
final.to_csv(output_file, index=False)
print(f"✅ Reformatted manifest with {len(final)} Tumor-Normal pairs saved to: {output_file}")

# Log unmatched cases
if unmatched:
    print("⚠️ Skipped the following case_id(s) due to missing Tumor or Normal:")
    for case_id in unmatched:
        print(f" - {case_id}")

