#!/bin/bash
#SBATCH --job-name=dedup_cram
#SBATCH --time=72:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --output=samtools_dedup_%j.out
#SBATCH --error=samtools_dedup_%j.err

set -euo pipefail
module load samtools

cd /data/Sherlock_Lung/CalebHartman/SV_calling/cram-no_duplicates

# Define input CRAMs
TUMOR_INPUT="/data/ITEB_Lung_WGS/FireCloud_downloads_2024/NSLC-AQFL-TTP1-A-1-1-D-A97M-36/NSLC-AQFL-TTP1-A-1-1-D-A97M-36.cram"
NORMAL_INPUT="/data/ITEB_Lung_WGS/FireCloud_downloads_2024/NSL01493_2000/NSL01493_2000.cram"

# Output BAMs
TUMOR_OUTPUT="NSLC-AQFL-TTP1-A-1-1-D-A97M-36.no_duplicates.bam"
NORMAL_OUTPUT="NSL01493_2000.no_duplicates.bam"

process_cram() {
    input_cram=$1
    output_bam=$2
    prefix=$(basename "${output_bam}" .no_duplicates.bam)

    echo "Processing: ${input_cram}"

    # Step 1: name sort (required for fixmate)
    samtools collate -@ 8 -o "${prefix}.namesort.bam" "${input_cram}"

    # Step 2: fixmate (adds ms/MC tags for markdup)
    samtools fixmate -m -@ 8 "${prefix}.namesort.bam" "${prefix}.fixmate.bam"

    # Step 3: coordinate sort for markdup
    samtools sort -@ 8 -T "${prefix}_tmp" -o "${prefix}.fixmate.sorted.bam" "${prefix}.fixmate.bam"

    # Step 4: mark duplicates
    samtools markdup -r -@ 8 "${prefix}.fixmate.sorted.bam" "${output_bam}"

    # Step 5: index final BAM
    samtools index "${output_bam}"

    # Step 6: cleanup
    rm -f "${prefix}.namesort.bam" "${prefix}.fixmate.bam" "${prefix}.fixmate.sorted.bam"
}

# Run for tumor and normal
process_cram "$TUMOR_INPUT" "$TUMOR_OUTPUT"
process_cram "$NORMAL_INPUT" "$NORMAL_OUTPUT"

echo "All steps completed successfully."
