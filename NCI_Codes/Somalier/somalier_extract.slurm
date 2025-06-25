#!/bin/bash
#SBATCH --job-name=somalier_extract
#SBATCH --output=/data/Sherlock_Lung/CalebHartman/420_Lung_Samples_Somalier/extracted_logs_2024/somalier_extract_%A_%a.out
#SBATCH --error=/data/Sherlock_Lung/CalebHartman/420_Lung_Samples_Somalier/extracted_logs_2024/somalier_extract_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --mem=16GB
#SBATCH --cpus-per-task=8
#SBATCH --array=1-410

module load somalier

# Set variables
WORK_DIR="/data/Sherlock_Lung/CalebHartman/420_Lung_Samples_Somalier"
SITES_VCF="${WORK_DIR}/sites.hg38.vcf.gz"
CRAM_LIST="${WORK_DIR}/cram_list_2024.txt"
OUTPUT_DIR="${WORK_DIR}/extracted_2024"
REFERENCE="/data/Sherlock_Lung/CalebHartman/Phuc_Conpair/Homo_sapiens_assembly38.fasta"

mkdir -p ${WORK_DIR}/extracted_logs_2024
mkdir -p ${OUTPUT_DIR}

# Read CRAM files and extract
CRAM_FILES=($(cat ${CRAM_LIST}))
CRAM_FILE="${CRAM_FILES[$((SLURM_ARRAY_TASK_ID - 1))]}"

somalier extract -d ${OUTPUT_DIR} --sites ${SITES_VCF} -f ${REFERENCE} ${CRAM_FILE}
