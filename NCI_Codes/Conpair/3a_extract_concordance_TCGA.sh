#!/bin/bash
#SBATCH --job-name=extract_concordance
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00


module load R
Rscript /data/Sherlock_Lung/CalebHartman/Phuc_Conpair/New_TCGA_Batch_2024/3a_extract_concordance_TCGA.R
