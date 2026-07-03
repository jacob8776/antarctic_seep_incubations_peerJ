#!/bin/bash

# Define Absolute Paths
FASTA_DIR="/nfs5/MICRO/Thurber_Lab/jacob/metagenome_experiments/antarctica_experiments/fastp_output/pear_output/fasta_files"
OUT_DIR="/nfs5/MICRO/Thurber_Lab/jacob/metagenome_experiments/antarctica_experiments/fastp_output/pear_output/prodigal_results"

mkdir -p "$OUT_DIR"

# Clear old task file
> prodigal_tasks.txt

echo "Generating compressed output commands..."

for f in "${FASTA_DIR}"/*.fasta; do
    name=$(basename "$f" .fasta)
    
    # We use a subshell ( ... ) to run prodigal and then immediately gzip the results.
    # 1. -f gff: Uses the much smaller GFF format instead of GBK.
    # 2. gzip: Compresses the FAA and FNA files as soon as they are finished.
    echo "prodigal -i $f -a ${OUT_DIR}/${name}.faa -d ${OUT_DIR}/${name}.fna -f gff -o ${OUT_DIR}/${name}.gff -p meta && gzip ${OUT_DIR}/${name}.faa ${OUT_DIR}/${name}.fna ${OUT_DIR}/${name}.gff" >> prodigal_tasks.txt
done

echo "Done! Created $(wc -l < prodigal_tasks.txt) commands with auto-compression."
