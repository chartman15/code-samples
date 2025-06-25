#!/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=16GB
#SBATCH --time=72:00:00

# Load required modules
module load manta
module load strelka

# Input arguments
tumor_bam=$1
normal_bam=$2
output_folder=$3

# Define temporary scratch space on local node
scratch="/lscratch/$SLURM_JOB_ID"

# Set up output and scratch directories
final_output="/data/Sherlock_Lung/CalebHartman/SV_calling/${output_folder}"
manta_dir="${scratch}/${output_folder}_manta"
strelka_dir="${scratch}/${output_folder}_strelka"

mkdir -p "${final_output}"
mkdir -p "${manta_dir}" "${strelka_dir}"

echo "Running Manta..."
configManta.py \
    --tumorBam "${tumor_bam}" \
    --normalBam "${normal_bam}" \
    --referenceFasta /data/Sherlock_Lung/Share/Reference/hg38/Homo_sapiens_assembly38.fasta \
    --runDir "${manta_dir}"

"${manta_dir}/runWorkflow.py" -m local -j "${SLURM_CPUS_PER_TASK}" -g $((SLURM_MEM_PER_NODE / 1024))

# Check for successful Manta run
manta_indels="${manta_dir}/results/variants/candidateSmallIndels.vcf.gz"
if [ ! -f "${manta_indels}" ]; then
    echo "ERROR: Manta failed to generate candidateSmallIndels.vcf.gz for ${output_folder}"
    exit 1
fi

cp -r "${manta_dir}/results" "${final_output}/manta_results"

echo "Running Strelka..."
configureStrelkaSomaticWorkflow.py \
    --tumorBam "${tumor_bam}" \
    --normalBam "${normal_bam}" \
    --callRegions /data/Sherlock_Lung/Share/Reference/hg38/callable.bed.gz \
    --indelCandidates "${manta_indels}" \
    --referenceFasta /data/Sherlock_Lung/Share/Reference/hg38/Homo_sapiens_assembly38.fasta \
    --runDir "${strelka_dir}"

"${strelka_dir}/runWorkflow.py" -m local -j "${SLURM_CPUS_PER_TASK}"

cp -r "${strelka_dir}/results" "${final_output}/strelka_results"

# Optionally clean up scratch space
# rm -rf "${manta_dir}" "${strelka_dir}"

