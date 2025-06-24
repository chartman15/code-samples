#!/bin/bash
set -euo pipefail

# Load required modules
module load angsd
module load R

# Input: path to BAM file
bam="$1"
prefix=$(basename "$bam" .bam)
outdir="$PWD/New_TCGA_2024_Admixture_fastNGSadmix/${prefix}"
#outdir="$PWD/New_Sherlock_Missing_Admixture_fastNGSadmix/${prefix}"
mkdir -p "$outdir"


# Use lscratch for fast I/O
scratchdir="/lscratch/${SLURM_JOB_ID}/${prefix}"
mkdir -p "$scratchdir"

# Reference and panel data
fastAD_install=/data/Sherlock_Lung/CalebHartman/Admixture/fastNGSadmix
SITES=${fastAD_install}/data1000genomes/data1000genomes/1000genomesRefPanel.sites
REF=${fastAD_install}/data1000genomes/data1000genomes/refPanel_1000genomesRefPanel.txt
NIND=${fastAD_install}/data1000genomes/data1000genomes/nInd_1000genomesRefPanel.txt
GENO=${fastAD_install}/data1000genomes/data1000genomes/1000genomesRefPanel

# Set threads
threads="${SLURM_CPUS_PER_TASK:-1}"

# Copy BAM locally to scratch
cp "$bam" "$scratchdir/"
cd "$scratchdir"

# Run ANGSD
angsd -i "$prefix.bam" -GL 2 -sites "$SITES" -doGlf 2 -doMajorMinor 3 -minMapQ 30 -minQ 20 \
      -doDepth 1 -doCounts 1 -out "$prefix.angsd_out" -P "$threads"

# Modify Beagle output
zcat "$prefix.angsd_out.beagle.gz" | sed 's/chr//g' > "$prefix.angsd_out.beagle2"

# Run fastNGSadmix
"$fastAD_install"/fastNGSadmix -likes "$prefix.angsd_out.beagle2" -fname "$REF" -Nname "$NIND" \
    -out "$prefix" -whichPops all

# Run PCA
Rscript "$fastAD_install"/R/fastNGSadmixPCA.R -likes "$prefix.angsd_out.beagle2" \
    -qopt "$prefix.qopt" -out "$prefix" -ref "$GENO"

# Move results back to output dir
mv "$prefix"* "$outdir/"

# Done
echo "Finished processing $bam"
