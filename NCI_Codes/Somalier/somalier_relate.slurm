#!/bin/bash
#SBATCH --job-name=somalier_relate
#SBATCH --output=/data/Sherlock_Lung/CalebHartman/420_Lung_Samples_Somalier/extracted_logs_2024/somalier_relate.out
#SBATCH --error=/data/Sherlock_Lung/CalebHartman/420_Lung_Samples_Somalier/extracted_logs_2024/somalier_relate.err
#SBATCH --time=04:00:00
#SBATCH --mem=32GB
#SBATCH --cpus-per-task=4

module load somalier

# Set variables
WORK_DIR="/data/Sherlock_Lung/CalebHartman/420_Lung_Samples_Somalier/Batch_2024"
OUTPUT_DIR="${WORK_DIR}/extracted_2024"
GROUPS_FILE="${WORK_DIR}/tumor_normal_groups_somalier.txt"
#GROUPS_FILE="${WORK_DIR}/tumor_normal_groups.txt"
#GROUPS_FILE="${WORK_DIR}/sex_groups.txt"
#GROUPS_FILE="${WORK_DIR}/female_groups.txt"
#GROUPS_FILE="${WORK_DIR}/male_groups.txt"
PED_FILE="${WORK_DIR}/somalier_manifest_mapped.ped"


# Wait until all .somalier files are ready (timeout: ~5 min)
for i in {1..30}; do
  COUNT=$(ls ${OUTPUT_DIR}/*.somalier 2>/dev/null | wc -l)
  if [[ $COUNT -eq 410 ]]; then
    echo "All 410 .somalier files found. Proceeding to relate."
    somalier relate --groups ${GROUPS_FILE} --ped ${PED_FILE} ${OUTPUT_DIR}/*.somalier
    #somalier relate ${OUTPUT_DIR}/*.somalier
    exit 0
  else
    echo "Only found $COUNT files. Waiting 10s... (attempt $i)"
    sleep 10
  fi
done

echo "Error: Not all .somalier files found after waiting. Found $COUNT."
exit 1
