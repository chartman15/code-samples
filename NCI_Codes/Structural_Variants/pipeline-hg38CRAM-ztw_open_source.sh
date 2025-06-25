#!/bin/bash
# DNAseq variant calling pipeline using open-source tools

set -euo pipefail
set -x

ml samtools; ml bazam; ml picard; ml gatk/4.4.0.0

bamfile=$1
foldname=$2
nt=$3
shortReadremove=${4:-false}

barcode=$(basename "$bamfile")
barcode=${barcode%%.bam}
barcode=${barcode%%.cram}

ulimit -s 10240

workdir="/lscratch/$SLURM_JOB_ID/$foldname"
workdir2="$(pwd)/$foldname"
mkdir -p "$workdir" "$workdir2"
logfile="$workdir2/run.log"
exec >"$logfile" 2>&1
cd "$workdir"

rush_path=/home/hartmancas/bin/rush

# Extract RG info
samtools view -H "$bamfile" | grep "^@RG" | \
  sed -e 's/\tDS:[^\t]*\([\t\.]\)/\t/' -e 's/\t\t*/\t/g' -e 's/\t$//' | \
  sed 's/\t/\\t/g' | sed 's/ /_/g' > RG_info.txt

# Split and convert
samtools split -@ "$nt" -u norg_file.bam "$bamfile"

ls ${barcode}_*.bam | "$rush_path" -k -v nt="$nt" \
  'samtools index -@ {nt} {}'

ls ${barcode}_*.bam | "$rush_path" -k -v nt="$nt" \
  'java -jar $BAZAMPATH/bazam.jar -n {nt} -bam {} | bgzip >{.}.fq.gz && rm -f {} {}.bai'

rm -f ${barcode}_*.bam*

if [ "$shortReadremove" = true ]; then
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

fasta="/data/Sherlock_Lung/Share/Reference/hg38/Homo_sapiens_assembly38.fasta"
dbsnp="/fdb/GATK_resource_bundle/hg38/dbsnp_138.hg38.vcf.gz"
known_Mills_indels="/fdb/GATK_resource_bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
known_1000G_indels="/fdb/GATK_resource_bundle/hg38/Homo_sapiens_assembly38.known_indels.vcf.gz"

# BWA Mapping + Sorting
bam_input_list=""
itot=$(ls ${barcode}_*.fastq.gz | grep -v "lowq.fastq.gz" | wc -l)
itot=$((itot - 1))

for i in $(seq 0 "$itot"); do
  fq="${barcode}_${i}.fastq.gz"
  rg=$(awk -v line="$i" 'NR==(line+1)' RG_info.txt)

  out_bam="sorted_set${i}.bam"
  bwa mem -M -t "$nt" -R "$rg" "$fasta" "$fq" | \
    samtools sort -@ "$nt" -o "$out_bam"
  bam_input_list="$bam_input_list I=$out_bam"
  rm -f "$fq"
done

# Mark duplicates
picard MarkDuplicates $bam_input_list \
  O=deduped.bam M=dedup_metrics.txt CREATE_INDEX=true

rm -f sorted_set*.bam*

# Coverage Metrics
picard CollectWgsMetrics \
  I=deduped.bam O=coverage_metrics.txt R="$fasta"

# Base Quality Recalibration
gatk BaseRecalibrator \
  -I deduped.bam -R "$fasta" --known-sites "$dbsnp" \
  --known-sites "$known_Mills_indels" --known-sites "$known_1000G_indels" \
  -O recal_data.table

gatk ApplyBQSR \
  -R "$fasta" -I deduped.bam --bqsr-recal-file recal_data.table \
  -O "$workdir2/${barcode}.cram" --create-output-bam-index true

rm -f deduped.bam* recal_data.table

# Move remaining outputs
mv RG_info.txt norg_file.bam "$workdir2/"
[ "$shortReadremove" = true ] && mv *_lowq.fastq.gz "$workdir2/"

ls -lh
