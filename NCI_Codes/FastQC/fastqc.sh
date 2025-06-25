#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH --mem=8g
#SBATCH --cpus-per-task=4
#SBATCH --output=fastqc_%j.out
#SBATCH --error=fastqc_%j.err

module load fastqc
mkdir -p fastqc_ERR194160_results
echo "Starting FastQC at $(date)"

fastqc ERR194160_NA12891_C0JVFACXX_*.fastq.gz -o fastqc_ERR194160_results/

echo "FastQC completed at $(date)"
echo "FastQC job finished successfully!" | tee -a fastqc_%j.out
