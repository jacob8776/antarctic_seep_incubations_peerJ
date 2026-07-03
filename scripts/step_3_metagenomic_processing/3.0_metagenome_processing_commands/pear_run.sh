#!/bin/bash

# Input and output directories
input_dir="/nfs5/MICRO/Thurber_Lab/jacob/metagenome_experiments/antarctica_experiments/fastp_output/"
output_dir="/nfs5/MICRO/Thurber_Lab/jacob/metagenome_experiments/antarctica_experiments/fastp_output/pear_output/"

# Make sure the output directory exists
mkdir -p "$output_dir"

# Loop through each forward read file
for forward_read in ${input_dir}*R1.fastq; do
    # Construct the reverse read file name by replacing R1 with R2
    reverse_read=${forward_read/R1.fastq/R2.fastq}

    # Extract the base name for output prefix (removing the directory path and file extension)
    base_name=$(basename "$forward_read" _R1.fastq)

    # Run pear command
    echo "pear -f "$forward_read" -r "$reverse_read" -p 0.001 -o "${output_dir}${base_name}" -y 64G -j 16"
done

