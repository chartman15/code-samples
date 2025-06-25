#! /bin/bash
module load cellranger || exit 1
## uncomment the following line if encountering 'resource unavailable' errors
## despite using --localcores and --localmem
# ulimit -u 4096
cellranger count --id=C4_SI_TT_G3 \
		 --transcriptome=/data/Sherlock_Lung/CalebHartman/cellranger_mkfastq/refdata-gex-GRCh38-2020-A \
		 --fastqs=/data/Sherlock_Lung/CalebHartman/cellranger_mkfastq/sample_C4_bcl_to_fastq/outs/fastq_path \
		 --sample=C4 \
		 --create-bam=false \
		 --localcores=8 \
		 --localmem=64

