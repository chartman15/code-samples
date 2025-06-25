#!/bin/bash

# Define output directory for FastQC results
OUTPUT_DIR="fastqc_results"

# Make sure the output directory exists
mkdir -p $OUTPUT_DIR

# Generate the swarm file
SWARM_FILE="fastqc.swarm"
> $SWARM_FILE  # Clear previous swarm file

# Loop over each FastQ file and add a command to the swarm file
for file in ERR194160_NA12891_C0JVFACXX_*_*.fastq.gz
do
    # Extract the directory (if needed) and filename
    DIR=$(dirname "$file")
    BASENAME=$(basename "$file")

    # Append the command to the swarm file with module load
    echo "module load fastqc; cd $DIR; fastqc -o $OUTPUT_DIR $BASENAME" >> $SWARM_FILE
done

echo "Swarm file created: $SWARM_FILE"
