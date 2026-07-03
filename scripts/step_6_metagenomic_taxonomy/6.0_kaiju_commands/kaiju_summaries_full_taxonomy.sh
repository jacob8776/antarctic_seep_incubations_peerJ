#!/bin/bash

# --- ⚙️ Configuration ---
# Define the input and output directories relative to where you run the script.
INPUT_DIR="/home/micro/wynneja/jacob/metagenome_experiments/antarctica_experiments/kaiju_annotation/short_annotations"  # Change this to your input folder path (e.g., "../raw_data/kaiju")
OUTPUT_DIR="/home/micro/wynneja/jacob/metagenome_experiments/antarctica_experiments/kaiju_annotation/short_annotations/genus_full_taxonomy" # Change this to your desired output folder path

# Define the common database files path
NODES_DMP="../kaiju_annotation/nodes.dmp"
NAMES_DMP="../kaiju_annotation/names.dmp"
RANK="genus" # The taxonomic rank to summarize

# --- 🚀 Setup and Processing ---

# 1. Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# 2. Loop over all files ending in .out in the specified input directory
# Note: We use "$INPUT_DIR"/*.out to ensure the search is constrained.
for INFILE_PATH in "$INPUT_DIR"/*.out
do
    # Check if the pattern matched any files (if no files, $INFILE_PATH remains "$INPUT_DIR"/*.out)
    if [ -f "$INFILE_PATH" ]; then
        
        # Extract just the filename from the full path (e.g., lane1-s001...kaiju.out)
        INFILE_NAME=$(basename "$INFILE_PATH")
        
        # Create the desired output filename by removing the ".out" extension 
        # and adding the "_summary.tsv" suffix (e.g., lane1-s001...kaiju_summary.tsv)
        BASE_NAME="${INFILE_NAME%.out}"
        OUTFILE_NAME="${BASE_NAME}_summary.tsv"
        
        # Define the full path for the output file
        OUTFILE_PATH="$OUTPUT_DIR/$OUTFILE_NAME"
        
        # Print the command (for review)
        echo "kaiju2table -t $NODES_DMP -n $NAMES_DMP -r $RANK -l superkingdom,phylum,class,order,family,genus -o $OUTFILE_PATH $INFILE_PATH"
        
    fi
done
