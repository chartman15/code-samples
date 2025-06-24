#!/bin/bash
set -euo pipefail

# Load modules
module load samtools bwa GATK

# Inputs
cram=$1
basename=$(basename "$cram" .cram)
#outdir=/data/Sherlock_Lung/CalebHartman/realign_hg19_output/New_Sherlock_Missing/$basename
outdir=/data/Sherlock_Lung/CalebHartman/realign_hg19_output/New_TCGA_LUAD/$basename
mkdir -p $outdir

# Reference genome (must match panel — GRCh37/hg19)
#ref=/data/Sherlock_Lung/Share/Reference/hg19/human_g1k_v37.fasta
ref=/data/Sherlock_Lung/Share/Reference/hg19/human_g1k_v37_decoy.fasta

# Work in lscratch
tmpdir=/lscratch/$SLURM_JOB_ID
mkdir -p $tmpdir
cd $tmpdir

# Copy CRAM to lscratch
cp $cram .

# Convert to FASTQ
samtools fastq -@ 8 -1 ${basename}_R1.fastq.gz -2 ${basename}_R2.fastq.gz \
    -0 /dev/null -s /dev/null -n $(basename $cram)

# Align to hg19
bwa mem -t 8 -R "@RG\tID:$basename\tSM:$basename\tPL:ILLUMINA" \
    $ref ${basename}_R1.fastq.gz ${basename}_R2.fastq.gz |
    samtools sort -@ 4 -o ${basename}.sorted.bam

# Mark duplicates
gatk --java-options "-Xmx16g" MarkDuplicates \
    -I ${basename}.sorted.bam \
    -O ${basename}.dedup.bam \
    -M ${basename}.metrics.txt

# Index
samtools index ${basename}.dedup.bam

# Copy results back to output dir
cp ${basename}.dedup.bam* ${basename}.metrics.txt $outdir/

# Clean up happens automatically when job ends
