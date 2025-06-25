#!/bin/bash
#SBATCH --job-name=Conpair
#SBATCH --partition=norm
#SBATCH --cpus-per-task=12      
#SBATCH --mem=80G               
#SBATCH --time=96:00:00         
#SBATCH --array=1-1156            
#SBATCH --output=logs/conpair_%A_%a.out
#SBATCH --error=logs/conpair_%A_%a.err

export CONPAIR_DIR=/usr/local/apps/conpair/0.2
export GATK_JAR=/gatk/gatk-package-4.6.0.0-local.jar
export PYTHONPATH=/usr/local/apps/conpair/0.2/scripts:/usr/local/apps/conpair/0.2/modules

module load conpair/0.2 
module load python 

filename_manifest=/data/Sherlock_Lung/CalebHartman/Phuc_Conpair/New_TCGA_Batch_2024/formatted_tcga_manifest_2024.csv

# Extract the line corresponding to this array task ID (skip header)
line=$(tail -n +2 "$filename_manifest" | sed -n "${SLURM_ARRAY_TASK_ID}p")

# Extract columns from CSV line
case_id=$(echo "$line" | cut -d',' -f1)
sample=$(echo "$line" | cut -d',' -f2)
filename_bam_tumour=$(echo "$line" | cut -d',' -f3)
filename_bam_normal=$(echo "$line" | cut -d',' -f4)

# DEBUG output to verify variables:
echo "DEBUG:"
echo "SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
echo "line=${line}"
echo "case_id=${case_id}"
echo "sample=${sample}"
echo "filename_bam_tumour=${filename_bam_tumour}"
echo "filename_bam_normal=${filename_bam_normal}"

ref=/data/Sherlock_Lung/CalebHartman/Phuc_Conpair/Homo_sapiens_assembly38.fasta
markers=/data/Sherlock_Lung/CalebHartman/Phuc_Conpair/GRCh38.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.liftover.txt
dir_output=/data/Sherlock_Lung/CalebHartman/Phuc_Conpair/New_TCGA_Batch_2024/output/${sample}

# Check if tumor and normal BAM files exist before proceeding
if [[ -f "$filename_bam_tumour" && -f "$filename_bam_normal" ]]; then
    echo "Tumor and Normal BAM files found. Proceeding..."

    # Create output directory if it doesn't exist
    if [ ! -d "$dir_output" ]; then
        mkdir -p "$dir_output"
    fi

    # Clean pileup files (only if they exist to avoid errors)
    if [[ -f "${dir_output}/${sample}_tumor_pileup" ]]; then
        awk -v FS=" " 'NF>=7' "${dir_output}/${sample}_tumor_pileup" > "${dir_output}/tmp.csv" && mv -f "${dir_output}/tmp.csv" "${dir_output}/${sample}_tumor_pileup"
    else
        echo "Warning: Tumor pileup file not found: ${dir_output}/${sample}_tumor_pileup"
    fi

    if [[ -f "${dir_output}/${sample}_normal_pileup" ]]; then
        awk -v FS=" " 'NF>=7' "${dir_output}/${sample}_normal_pileup" > "${dir_output}/tmp.csv" && mv -f "${dir_output}/tmp.csv" "${dir_output}/${sample}_normal_pileup"
    else
        echo "Warning: Normal pileup file not found: ${dir_output}/${sample}_normal_pileup"
    fi

    # Run Conpair concordance and contamination estimation
    ${CONPAIR_DIR}/scripts/verify_concordance.py -T "${dir_output}/${sample}_tumor_pileup" -N "${dir_output}/${sample}_normal_pileup" -M "$markers" -O "${dir_output}/${sample}_concordance.txt"
    ${CONPAIR_DIR}/scripts/verify_concordance.py -T "${dir_output}/${sample}_tumor_pileup" -N "${dir_output}/${sample}_normal_pileup" -M "$markers" -H -O "${dir_output}/${sample}_homozygous_concordance.txt"
    ${CONPAIR_DIR}/scripts/estimate_tumor_normal_contamination.py -T "${dir_output}/${sample}_tumor_pileup" -N "${dir_output}/${sample}_normal_pileup" -M "$markers" -O "${dir_output}/${sample}_cross_contamination.txt"

else
    echo "Tumor or Normal BAM file not found for sample $sample."
    echo "Tumor BAM: $filename_bam_tumour"
    echo "Normal BAM: $filename_bam_normal"
fi
