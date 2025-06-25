#!/bin/bash
# swarm_script.sh

# Create a log file for tracking errors
LOG_FILE="fastp_trim.log"

# Create an empty swarm file
swarm_file="fastp.swarm"
> $swarm_file  # Clear the file if it exists

# Loop over all files matching the naming pattern for paired-end reads
for fastq1 in *_1.fastq.gz; do
    # Get the corresponding _2 file for paired-end reads
    fastq2="${fastq1/_1.fastq.gz/_2.fastq.gz}"

    # Check if the _2 file exists for paired-end reads
    if [[ -f "$fastq2" ]]; then
        # Output file names
        output_prefix="${fastq1/_1.fastq.gz/}"
        output_file1="${output_prefix}_1.trimmed.fastq.gz"
        output_file2="${output_prefix}_2.trimmed.fastq.gz"
        json_output="${output_prefix}.json"

        # Write the fastp command to the swarm file, including module loading
        echo "module load fastp; fastp -i $fastq1 -I $fastq2 -o $output_file1 -O $output_file2 -j $json_output -h ${output_prefix}.html" >> $swarm_file
    else
        # If the _2 file doesn't exist, print an error message to the log file
        echo "ERROR: Pair for $fastq1 not found" >> $LOG_FILE
    fi
done

echo "Swarm file generated: $swarm_file"
echo "To submit the job, run: swarm -f $swarm_file -g 8 -t 1 --logdir ./swarm_logs"
