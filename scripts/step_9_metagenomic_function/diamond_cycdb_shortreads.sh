# 1. Define the directories
INPUT_DIR="/nfs5/MICRO/Thurber_Lab/jacob/metagenome_experiments/antarctica_experiments/fastp_output/pear_output/prodigal_results"
REF_DB="/nfs5/MICRO/Thurber_Lab/jacob/metagenome_experiments/antarctica_experiments/diamond_cycdb/combined_cyc.dmnd"
OUT_DIR="/nfs5/MICRO/Thurber_Lab/jacob/metagenome_experiments/antarctica_experiments/diamond_cycdb/short_read_analysis"

# 2. Make sure the output folder exists
mkdir -p "$OUT_DIR"

# 3. Clear the task file
> diamond_tasks_short_reads.txt

# 4. The Loop that creates 34 specific lines
for faa in "${INPUT_DIR}"/*.faa; do
    # This part extracts the actual sample name (e.g., lane1-s001...)
    sample=$(basename "$faa" .faa)
    
    # This writes a specific line for EACH sample into the text file
    echo "diamond blastp --threads 16 --db ${REF_DB} --query ${faa} --out ${OUT_DIR}/${sample}_vs_cycdb.tsv --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle --evalue 1e-5 --max-target-seqs 1 --block-size 12.0 --sensitive --id 50 --query-cover 50" >> diamond_tasks_short_reads.txt
done

# 5. Check the result
echo "Created $(wc -l < diamond_tasks_short_reads.txt) specific commands."
head -n 1 diamond_tasks_short_reads.txt
