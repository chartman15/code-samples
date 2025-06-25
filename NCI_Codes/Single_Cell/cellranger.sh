#! /bin/bash
module load cellranger || exit 1
## uncomment the following line if encountering 'resource unavailable' errors 
## despite using --localcores and --localmem
# ulimit -u 4096
cellranger mkfastq --id=sample_C4_bcl_to_fastq \
		   --run=/data/Sherlock_Lung/scRNAseq/Chongyi_pilot/240415_VH01509_7_AAFLMLVM5/ \
        	   --csv=/data/Sherlock_Lung/CalebHartman/cellranger_mkfastq/Sample_Sheet_C4_SI_TT_G3.csv \
        	   --localcores=$SLURM_CPUS_PER_TASK \
        	   --localmem=34

