#!/bin/bash
# *******************************************
# Script to perform DNA seq variant calling
# using a single sample with multiple input FASTQ sets
# *******************************************
set -euo pipefail
set -x  # Enable command tracing for debugging

ml samtools; ml bazam

# Input variables (from GitHub usage instructions)
bamfile=$1           # sample.cram
foldname=$2          # SampleOutput
nt=$3                # 16
shortReadremove=${4:-false}  # true/false

barcode=$(basename "$bamfile")
barcode=${barcode%%.bam}
barcode=${barcode%%.cram}

ulimit -s 10240

# ******************************************
# 0. Setup
# ******************************************
echo "Setting up work directories..."
workdir="/lscratch/$SLURM_JOB_ID/$foldname"
workdir2="$(pwd)/$foldname"
mkdir -p "$workdir" "$workdir2"
logfile="$workdir2/run.log"
exec >"$logfile" 2>&1
cd "$workdir"

rush_path=/home/hartmancas/bin/rush

# Extract clean RG info
echo "Extracting read group info from BAM file..."
samtools view -H "$bamfile" | grep "^@RG" | \
  sed -e 's/\tDS:[^\t]*\([\t\.]\)/\t/' -e 's/\t\t*/\t/g' -e 's/\t$//' | \
  sed 's/\t/\\t/g' | sed 's/ /_/g' > RG_info.txt

# Split BAM by read group
echo "Splitting BAM file by read group..."
#samtools split -@ "$nt" "$bamfile"
samtools split -@ "$nt" -u norg_file.bam "$bamfile"

# Index BAMs
echo "Indexing BAM files..."
ls ${barcode}_*.bam | "$rush_path" --dry-run -k -v nt="$nt" 'samtools index -@ {nt} {}' | sh

# Convert BAMs to FASTQ
echo "Converting BAM files to FASTQ..."
ls ${barcode}_*.bam | "$rush_path" --dry-run -k -v nt="$nt" \
  'java -jar $BAZAMPATH/bazam.jar -n {nt} -bam {} | bgzip >{.}.fq.gz && rm -f {} {}.bai' | sh

rm -f ${barcode}_*.bam*

# Optional: Remove short reads < 30bp
if [ "$shortReadremove" = true ]; then
  echo "Removing short reads < 30bp..."
  ls *.fq.gz | "$rush_path" -k -v p="'" \
    'zcat {} | paste - - - - - - - - | \
    awk -F "\t" -v OFS="\t" '"$p"'{
      if(length($2) > 30 && length($6) > 30)
        print $1"\n"$2"\n"$3"\n"$4"\n"$5"\n"$6"\n"$7"\n"$8 > "/dev/stdout";
    }'"$p"' | bgzip >{..}_lowq.fastq.gz'
  rm -f *.fq.gz
else
  for file in *.fq.gz; do
    mv "$file" "${file%.fq.gz}.fastq.gz"
  done
fi

ls -hl

# Reference and Sentieon paths
fasta="/data/Sherlock_Lung/Share/Reference/hg38/Homo_sapiens_assembly38.fasta"
dbsnp="/fdb/GATK_resource_bundle/hg38/dbsnp_138.hg38.vcf.gz"
known_Mills_indels="/fdb/GATK_resource_bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
known_1000G_indels="/fdb/GATK_resource_bundle/hg38/Homo_sapiens_assembly38.known_indels.vcf.gz"
release_dir="/data/hartmancas/sentieon-genomics-202308.03"

export SENTIEON_INSTALL_DIR=$release_dir
export SENTIEON_LICENSE=badmin.cit.nih.gov:18353
export SENTIEON_AUTH_DATA="26Mu6qUhnH+HjRptf0EZ6ii7Aqw="

# ******************************************
# 1. Mapping FASTQ sets with BWA-MEM, sorting
# ******************************************
echo "Starting BWA-MEM mapping and sorting..."
bam_input=""
itot=$(ls ${barcode}_*.fastq.gz | grep -v "lowq.fastq.gz" | wc -l)
itot=$((itot - 1))

for i in $(seq 0 "$itot"); do
  filei="${barcode}_${i}.fastq.gz"
  rgi=$(awk -v line="$i" 'NR==(line+1)' RG_info.txt)
  bam_input="$bam_input -i sorted_set${i}.bam"

  echo "Processing $filei with RG: $rgi"

  ( "$SENTIEON_INSTALL_DIR/bin/sentieon" bwa mem -M -p -R "$rgi" -t "$nt" -K 10000000 "$fasta" "$filei" || echo "BWA-MEM error for $filei" ) \
    | "$SENTIEON_INSTALL_DIR/bin/sentieon" util sort -r "$fasta" -o "sorted_set${i}.bam" -t "$nt" --sam2bam -i -

  rm -f "$filei"
done

ls -hl

# ******************************************
# 3. Remove duplicates
# ******************************************
echo "Removing duplicates..."
"$SENTIEON_INSTALL_DIR/bin/sentieon" driver -t "$nt" $bam_input --algo LocusCollector --fun score_info score.txt
"$SENTIEON_INSTALL_DIR/bin/sentieon" driver -t "$nt" $bam_input --algo Dedup --score_info score.txt --metrics dedup_metrics.txt "$workdir/deduped.bam"
#"$SENTIEON_INSTALL_DIR/bin/sentieon" driver -t "$nt" $bam_input --algo Dedup --rmdup --score_info score.txt --metrics dedup_metrics.txt "$workdir/deduped.bam"

rm -f sorted_set*.bam*

# ******************************************
# 4. Coverage metrics
# ******************************************
echo "Calculating coverage metrics..."
"$SENTIEON_INSTALL_DIR/bin/sentieon" driver -r "$fasta" -t "$nt" -i "$workdir/deduped.bam" --algo CoverageMetrics coverage_metrics

# ******************************************
# 5. Realignment
# ******************************************
echo "Realigning BAM file..."
"$SENTIEON_INSTALL_DIR/bin/sentieon" driver -r "$fasta" -t "$nt" -i "$workdir/deduped.bam" --algo Realigner \
  -k "$known_Mills_indels" -k "$known_1000G_indels" "$workdir/realigned.bam"
rm -f "$workdir/deduped.bam"

# ******************************************
# 6. Base Recalibration
# ******************************************
echo "Performing base recalibration..."
"$SENTIEON_INSTALL_DIR/bin/sentieon" driver -r "$fasta" -t "$nt" -i "$workdir/realigned.bam" --algo QualCal \
  -k "$dbsnp" -k "$known_Mills_indels" -k "$known_1000G_indels" "$workdir/recal_data.table"

"$SENTIEON_INSTALL_DIR/bin/sentieon" driver -r "$fasta" -t "$nt" -i "$workdir/realigned.bam" \
  --read_filter QualCalFilter,table="$workdir/recal_data.table",indel=false,levels=10/20/30/40/50 \
  --algo ReadWriter "$workdir2/${barcode}.cram"

rm -f "$workdir/realigned.bam"*
rm -f "$workdir/recal_data.table"*

# ******************************************
# Final Steps
# ******************************************
echo "Moving final output files..."
#mv RG_info.txt "$workdir2/"
mv RG_info.txt norg_file.bam "$workdir2/"
[ "$shortReadremove" = true ] && mv *_lowq.fastq.gz "$workdir2/"

ls -lh
