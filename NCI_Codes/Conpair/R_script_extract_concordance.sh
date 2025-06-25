#!/bin/bash
#SBATCH --job-name=R_conpair_extract
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00

module load R
Rscript /data/Sherlock_Lung/CalebHartman/Phuc_Conpair/3a_extract_concordance_Sherlock.R
