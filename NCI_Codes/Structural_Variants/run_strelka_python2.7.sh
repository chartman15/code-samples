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

mkdir -p "${final_output}" "${manta_dir}" "${strelka_dir}"

echo "===== RUNNING MANTA WITH PYTHON 2.7 ====="

# Save original PATH and LD_LIBRARY_PATH
ORIG_PATH="$PATH"
ORIG_LD_LIBRARY_PATH="$LD_LIBRARY_PATH"

# Temporarily prepend Python 2.7 to PATH and LD_LIBRARY_PATH
export PATH="$HOME/bin/python2.7/install/bin:$PATH"
export LD_LIBRARY_PATH="$HOME/bin/python2.7/install/lib:$LD_LIBRARY_PATH"

# Run Manta
configManta.py \
    --tumorBam "${tumor_bam}" \
    --normalBam "${normal_bam}" \
    --referenceFasta /data/Sherlock_Lung/Share/Reference/hg38/Homo_sapiens_assembly38.fasta \
    --runDir "${manta_dir}"

"${manta_dir}/runWorkflow.py" -m local -j "${SLURM_CPUS_PER_TASK}" -g $((SLURM_MEM_PER_NODE / 1024))

# Check Manta output
manta_indels="${manta_dir}/results/variants/candidateSmallIndels.vcf.gz"
if [ ! -f "${manta_indels}" ]; then
    echo "ERROR: Manta failed to generate candidateSmallIndels.vcf.gz for ${output_folder}"
    exit 1
fi

cp -r "${manta_dir}/results" "${final_output}/manta_results"

# Restore environment
export PATH="$ORIG_PATH"
export LD_LIBRARY_PATH="$ORIG_LD_LIBRARY_PATH"

echo "===== RUNNING STRELKA WITH SYSTEM PYTHON ====="

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
