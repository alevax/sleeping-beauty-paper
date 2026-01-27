#!/bin/bash
#SBATCH --job-name=salmon_quant
#SBATCH --output=logs/salmon_quant_%A_%a.log
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --partition=20
#SBATCH --array=1-43

#
# slurm_salmon_quant.sh
#
# Quantify transcript abundance using Salmon
# Handles both single-end and paired-end samples
#
# Usage: sbatch scripts/slurm_salmon_quant.sh (from isoform_switching/)
#

set -e

# Salmon is available at /usr/local/bin/salmon (v1.8.0)

THREADS=${SLURM_CPUS_PER_TASK:-8}

# === Configuration ===
FASTQ_DIR="data/fastq"
SAMPLES_FILE="data/fastq/samples.txt"
INDEX_DIR="reference/salmon_index"
OUTPUT_DIR="results/salmon"

mkdir -p "$OUTPUT_DIR"
mkdir -p logs

# Get sample info for this array task
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLES_FILE")
SAMPLE=$(echo "$LINE" | cut -d',' -f1)
SEQ_TYPE=$(echo "$LINE" | cut -d',' -f2)

echo "=== Salmon Quantification ==="
echo "Date: $(date)"
echo "Node: $(hostname)"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Sample: $SAMPLE"
echo "Sequencing type: $SEQ_TYPE"
echo "Threads: $THREADS"
echo ""

# Output directory for this sample
SAMPLE_OUT="$OUTPUT_DIR/$SAMPLE"

if [ -d "$SAMPLE_OUT" ] && [ -f "$SAMPLE_OUT/quant.sf" ]; then
    echo "SKIP: $SAMPLE_OUT already exists"
    exit 0
fi

# Find FASTQ files
FASTQ_R1=$(ls ${FASTQ_DIR}/${SAMPLE}_*_R1_*.fastq.gz 2>/dev/null | head -1)

if [ -z "$FASTQ_R1" ]; then
    echo "ERROR: No FASTQ found for $SAMPLE"
    exit 1
fi

echo "R1: $FASTQ_R1"

# === Run Salmon ===
if [ "$SEQ_TYPE" == "paired" ]; then
    FASTQ_R2=$(ls ${FASTQ_DIR}/${SAMPLE}_*_R2_*.fastq.gz 2>/dev/null | head -1)

    if [ -z "$FASTQ_R2" ]; then
        echo "ERROR: Missing R2 for paired-end sample $SAMPLE"
        exit 1
    fi

    echo "R2: $FASTQ_R2"
    echo ""
    echo "Running Salmon (paired-end)..."

    salmon quant \
        -i "$INDEX_DIR" \
        -l A \
        -1 "$FASTQ_R1" \
        -2 "$FASTQ_R2" \
        -p $THREADS \
        --validateMappings \
        --gcBias \
        --seqBias \
        -o "$SAMPLE_OUT"
else
    echo ""
    echo "Running Salmon (single-end)..."

    salmon quant \
        -i "$INDEX_DIR" \
        -l A \
        -r "$FASTQ_R1" \
        -p $THREADS \
        --validateMappings \
        --gcBias \
        --seqBias \
        -o "$SAMPLE_OUT"
fi

echo ""
echo "Output: $SAMPLE_OUT"
echo "Date: $(date)"
echo "=== Done ==="
