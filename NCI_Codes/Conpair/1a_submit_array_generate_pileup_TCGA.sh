#!/bin/bash
#SBATCH --job-name=conpair_pileup
#SBATCH --partition=norm
#SBATCH --cpus-per-task=12      
#SBATCH --mem=80G               
#SBATCH --time=96:00:00         
#SBATCH --array=1-1156
#SBATCH --output=/data/Sherlock_Lung/CalebHartman/Phuc_Conpair/New_TCGA_Batch_2024/logs/conpair_pileup_%A_%a.out
#SBATCH --error=/data/Sherlock_Lung/CalebHartman/Phuc_Conpair/New_TCGA_Batch_2024/logs/conpair_pileup_%A_%a.err

module load GATK

filename_manifest=/data/Sherlock_Lung/CalebHartman/Phuc_Conpair/New_TCGA_Batch_2024/formatted_tcga_manifest_2024.csv

# Skip the header line using tail
line=$(tail -n +2 "$filename_manifest" | sed -n "${SLURM_ARRAY_TASK_ID}p")

# Extract columns with comma delimiter
case_id=$(echo "$line" | cut -d',' -f1)
sample=$(echo "$line" | cut -d',' -f2)
filename_bam_tumour=$(echo "$line" | cut -d',' -f3)
filename_bam_normal=$(echo "$line" | cut -d',' -f4)

ref=/data/Sherlock_Lung/CalebHartman/Phuc_Conpair/Homo_sapiens_assembly38.fasta
markers=/data/Sherlock_Lung/CalebHartman/Phuc_Conpair/GRCh38.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.liftover.bed
dir_output=/data/Sherlock_Lung/CalebHartman/Phuc_Conpair/New_TCGA_Batch_2024/output/${sample}

# Run only if both BAM files exist
if [[ -f $filename_bam_normal && -f $filename_bam_tumour ]]; then
  # Create output directory if it doesn't exist
  mkdir -p $dir_output
  
  # Run GATK Pileup for tumor and normal BAM files
  gatk Pileup -R $ref -L $markers -I $filename_bam_tumour -verbose true -O ${dir_output}/${sample}_tumor_pileup
  gatk Pileup -R $ref -L $markers -I $filename_bam_normal -verbose true -O ${dir_output}/${sample}_normal_pileup
else
  echo "Skipping sample $sample: Tumor or Normal BAM file not found."
fi
