#!/bin/bash

# 1. Base directories
input_base_dir="/nfs5/MICRO/Thurber_Lab/jacob/metagenome_experiments/antarctica_experiments/fastp_output"
# UPDATED OUTPUT DIRECTORY:
output_base_dir="/nfs5/MICRO/Thurber_Lab/jacob/metagenome_experiments/antarctica_experiments/kaiju_annotation/short_annotations"

# 2. Database paths
nodes="/nfs5/MICRO/Thurber_Lab/jacob/metagenome_experiments/antarctica_experiments/kaiju_annotation/nodes.dmp"
reference_db="/nfs5/MICRO/Thurber_Lab/jacob/metagenome_experiments/antarctica_experiments/kaiju_annotation/kaiju_db_refseq_nr.fmi"

# Create the new output directory if it doesn't exist
mkdir -p "$output_base_dir"

# 3. Loop through all R1 files
for r1_file in "$input_base_dir"/*_processed_R1.fastq; do
    
    # Check if file exists to prevent errors
    [ -e "$r1_file" ] || continue

    # Define the matching R2 file (Replace R1 with R2 in the filename)
    r2_file="${r1_file/_processed_R1.fastq/_processed_R2.fastq}"

    # Extract clean sample name for the output file
    basename=$(basename "$r1_file")
    sample_name="${basename%_processed_R1.fastq}"
    
    output_file="${output_base_dir}/${sample_name}_kaiju.out"

    # Print the Kaiju command
    echo "kaiju -z 16 -t ${nodes} -f ${reference_db} -i ${r1_file} -j ${r2_file} -o ${output_file}"

done
