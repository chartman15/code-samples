#!/bin/bash
#SBATCH --job-name=conpair_pileup
#SBATCH --partition=norm
#SBATCH --cpus-per-task=12      
#SBATCH --mem=80G               
#SBATCH --time=96:00:00         
#SBATCH --array=1-206            
#SBATCH --output=logs/conpair_pileup_%A_%a.out
#SBATCH --error=logs/conpair_pileup_%A_%a.err

module load GATK

filename_manifest=/data/Sherlock_Lung/CalebHartman/Phuc_Conpair/formatted_manifest_2024.csv
# Skip the header line using tail
line=$(tail -n +2 "$filename_manifest" | sed -n "${SLURM_ARRAY_TASK_ID}p")

# Extract columns
sample=$(echo "$line" | cut -d',' -f2)
filename_bam_tumour=$(echo "$line" | cut -d',' -f3)
filename_bam_normal=$(echo "$line" | cut -d',' -f4)

ref=/data/Sherlock_Lung/CalebHartman/Phuc_Conpair/Homo_sapiens_assembly38.fasta
markers=/data/Sherlock_Lung/CalebHartman/Phuc_Conpair/GRCh38.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.liftover.bed
dir_output=/data/Sherlock_Lung/CalebHartman/Phuc_Conpair/output/Sherlock/${sample}

#run
if [[ -f $filename_bam_normal && -f $filename_bam_tumour ]] ##Only run if match T/N files exist
then
  # Create all necessary directories if they do not exist
  if [ ! -d $dir_output ]; then
    mkdir -p $dir_output  # Create the full directory path
  fi
  
  # Run GATK Pileup command for tumor and normal BAM files
  gatk Pileup -R $ref -L $markers -I $filename_bam_tumour -verbose true -O ${dir_output}/${sample}_tumor_pileup
  gatk Pileup -R $ref -L $markers -I $filename_bam_normal -verbose true -O ${dir_output}/${sample}_normal_pileup
fi
