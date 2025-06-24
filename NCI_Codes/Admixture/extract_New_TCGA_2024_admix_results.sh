#!/bin/bash

base_dir="/data/Sherlock_Lung/CalebHartman/Admixture/New_TCGA_2024_Admixture_fastNGSadmix"
output="admixture_qopt_summary_with_ancestry.tsv"

# Population to superpopulation mapping
declare -A SUPERPOP=(
    [FIN]="EUR"
    [PEL]="AMR/Mixed"
    [PJL]="SAS"
    [CEU]="EUR"
    [YRI]="AFR"
    [CHB]="EAS"
)

printed_header=false
> "$output"

for qopt_file in "$base_dir"/*/*.qopt; do
    # Extract base filename without .qopt extension
    full_name=$(basename "$qopt_file" .qopt)

    # Try to extract a UUID
    uuid=$(echo "$full_name" | grep -oP '[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}')

    # If UUID not found, extract TCGA barcode (up to the first dot)
    if [[ -z "$uuid" ]]; then
        uuid="${full_name%%.*}"
    fi

    # Read header and values
    read -r header_line < "$qopt_file"
    read -r values_line < <(tail -n 1 "$qopt_file")

    # Split headers and values into arrays
    IFS=' ' read -ra headers <<< "$header_line"
    IFS=' ' read -ra values <<< "$values_line"

    # Find index of max value
    max_index=0
    max_value=${values[0]}
    for i in "${!values[@]}"; do
        if (( $(echo "${values[i]} > $max_value" | bc -l) )); then
            max_value=${values[i]}
            max_index=$i
        fi
    done

    # Find population and superpopulation
    max_pop=${headers[$max_index]}
    max_superpop=${SUPERPOP[$max_pop]}

    # Write header once
    if ! $printed_header; then
        echo -e "Sample_ID\t${header_line}\tMax_Ancestry" >> "$output"
        printed_header=true
    fi

    echo -e "$uuid\t${values_line}\t$max_superpop" >> "$output"
done
