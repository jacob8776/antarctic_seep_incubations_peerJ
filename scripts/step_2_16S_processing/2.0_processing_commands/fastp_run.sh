#!/bin/bash

# Directory containing your fastq.gz files
input_dir="/Users/wynne/Dropbox/masters_paper/updated_16s_run/unprocessed"
output_dir="/Users/wynne/Dropbox/masters_paper/updated_16s_run/raw"

# Make sure the output directory exists
mkdir -p "$output_dir"

# Loop over all the fastq.gz files in the input directory
for file in "$input_dir"/*.fastq.gz; do
  # Extract the base name of the file (without path and extension)
  base_name=$(basename "$file" .fastq.gz)
  
  # Construct the output file name
  output_file="$output_dir/${base_name}_fastp_output.fastq.gz"
  
  # Run fastp on the file
  fastp -i "$file" -o "$output_file"
done
