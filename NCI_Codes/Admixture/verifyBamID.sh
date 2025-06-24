#!/bin/bash

bam=$1
prefix=${bam%%.cram}
prefix=${prefix##*/}

ml verifybamid

verifybamid \
	--SVDPrefix /usr/local/apps/verifybamid/conda/envs/2.0.1/share/verifybamid2-2.0.1-8/resource/1000g.phase3.100k.b38.vcf.gz.dat \
	--BamFile ${bam} \
	--Reference /data/Sherlock_Lung/Share/Reference/hg38/Homo_sapiens_assembly38.fasta \
	--NumThread $SLURM_CPUS_PER_TASK --NumPC 4 --Output ${prefix}
