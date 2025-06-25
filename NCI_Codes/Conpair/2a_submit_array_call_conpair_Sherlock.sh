#!/bin/bash
#SBATCH --job-name=Conpair
#SBATCH --partition=norm
#SBATCH --cpus-per-task=12      
#SBATCH --mem=80G               
#SBATCH --time=96:00:00         
#SBATCH --array=1-206            
#SBATCH --output=logs/conpair_%A_%a.out
#SBATCH --error=logs/conpair_%A_%a.err

export CONPAIR_DIR=/usr/local/apps/conpair/0.2
export GATK_JAR=/gatk/gatk-package-4.6.0.0-local.jar
export PYTHONPATH=/usr/local/apps/conpair/0.2/scripts:/usr/local/apps/conpair/0.2/modules

module load conpair/0.2 
module load python 


filename_manifest=/data/Sherlock_Lung/CalebHartman/Phuc_Conpair/formatted_manifest_2024.csv
# Skip the header line using tail
line=$(tail -n +2 "$filename_manifest" | sed -n "${SLURM_ARRAY_TASK_ID}p")

# Extract columns
sample=$(echo "$line" | cut -d',' -f2)
filename_bam_tumour=$(echo "$line" | cut -d',' -f3)
filename_bam_normal=$(echo "$line" | cut -d',' -f4)

ref=/data/Sherlock_Lung/CalebHartman/Phuc_Conpair/Homo_sapiens_assembly38.fasta
markers=/data/Sherlock_Lung/CalebHartman/Phuc_Conpair/GRCh38.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.liftover.txt
dir_output=/data/Sherlock_Lung/CalebHartman/Phuc_Conpair/output/Sherlock/${sample}

#run
if [[ -f $filename_bam_normal && -f $filename_bam_tumour ]] ##Only run if match T/N files exist
  then
	if [ ! -d $dir_output ]
  		then
      		mkdir $dir_output
  	fi
  awk -v FS=" " "NF>=7" ${dir_output}/${sample}_tumor_pileup > ${dir_output}/tmp.csv && mv -f ${dir_output}/tmp.csv ${dir_output}/${sample}_tumor_pileup
  awk -v FS=" " "NF>=7" ${dir_output}/${sample}_normal_pileup > ${dir_output}/tmp.csv && mv -f ${dir_output}/tmp.csv ${dir_output}/${sample}_normal_pileup
  ${CONPAIR_DIR}/scripts/verify_concordance.py -T ${dir_output}/${sample}_tumor_pileup -N ${dir_output}/${sample}_normal_pileup -M $markers -O ${dir_output}/${sample}_concordance.txt
  ${CONPAIR_DIR}/scripts/verify_concordance.py -T ${dir_output}/${sample}_tumor_pileup -N ${dir_output}/${sample}_normal_pileup -M $markers -H -O ${dir_output}/${sample}_homozygous_concordance.txt 
  ${CONPAIR_DIR}/scripts/estimate_tumor_normal_contamination.py -T ${dir_output}/${sample}_tumor_pileup -N ${dir_output}/${sample}_normal_pileup -M $markers -O ${dir_output}/${sample}_cross_contamination.txt
fi
