#!/bin/bash
#SBATCH --job-name=align_chrSB
#SBATCH --output=logs/align_%A_%a.log
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=20
#SBATCH --array=1-29

#
# slurm_align_array.sh
#
# Align RNA-seq FASTQs to mm10+chrSB reference using HISAT2
# Uses SLURM array jobs to process all 29 paired-end samples in parallel
#
# Prerequisites:
#   - Run slurm_build_reference.sh first
#   - data/fastq/samples.txt file with sample info
#
# Usage: sbatch scripts/slurm_align_array.sh (from project root)
#

set -e  # Exit on error

THREADS=${SLURM_CPUS_PER_TASK:-8}

# === Setup ===
mkdir -p results/bams/stats

# === Load modules ===
module load samtools/1.17
module load hisat2/2.2.1

# Get sample info for this array task
# samples.txt format: SAMPLE_ID,TYPE (single or paired)
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" data/fastq/samples.txt)
SAMPLE=$(echo "$LINE" | cut -d',' -f1)
SEQ_TYPE=$(echo "$LINE" | cut -d',' -f2)

echo "=== Alignment Job Started ==="
echo "Date: $(date)"
echo "Node: $(hostname)"
echo "Working dir: $(pwd)"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Sample: $SAMPLE"
echo "Sequencing type: $SEQ_TYPE"
echo "Threads: $THREADS"
echo ""

# Output BAM
OUTBAM="results/bams/${SAMPLE}_sorted.bam"

if [ -f "$OUTBAM" ]; then
    echo "SKIP: $OUTBAM already exists"
    exit 0
fi

# === Align based on sequencing type ===
if [ "$SEQ_TYPE" == "single" ]; then
    # Single-end alignment
    FASTQ_R1=$(ls data/fastq/${SAMPLE}_*R1*.fastq.gz 2>/dev/null | head -1)

    if [ -z "$FASTQ_R1" ]; then
        echo "ERROR: No FASTQ found for $SAMPLE"
        exit 1
    fi

    echo "Input FASTQ: $FASTQ_R1"
    echo "Running HISAT2 (single-end)..."

    hisat2 -p $THREADS \
        -x reference/hisat2_index/mm10_chrSB \
        -U "$FASTQ_R1" \
        --dta \
        --summary-file "results/bams/stats/${SAMPLE}_alignment_summary.txt" \
        2> "results/bams/stats/${SAMPLE}_hisat2.log" \
    | samtools sort -@ $THREADS -o "$OUTBAM" -

else
    # Paired-end alignment
    FASTQ_R1=$(ls data/fastq/${SAMPLE}_*R1*.fastq.gz 2>/dev/null | head -1)
    FASTQ_R2=$(ls data/fastq/${SAMPLE}_*R2*.fastq.gz 2>/dev/null | head -1)

    if [ -z "$FASTQ_R1" ] || [ -z "$FASTQ_R2" ]; then
        echo "ERROR: Missing FASTQ pair for $SAMPLE"
        exit 1
    fi

    echo "Input R1: $FASTQ_R1"
    echo "Input R2: $FASTQ_R2"
    echo "Running HISAT2 (paired-end)..."

    hisat2 -p $THREADS \
        -x reference/hisat2_index/mm10_chrSB \
        -1 "$FASTQ_R1" \
        -2 "$FASTQ_R2" \
        --dta \
        --summary-file "results/bams/stats/${SAMPLE}_alignment_summary.txt" \
        2> "results/bams/stats/${SAMPLE}_hisat2.log" \
    | samtools sort -@ $THREADS -o "$OUTBAM" -
fi

# Index BAM
echo "Indexing BAM..."
samtools index "$OUTBAM"

# Report stats
echo ""
echo "=== Alignment Statistics ==="
TOTAL_READS=$(samtools view -c "$OUTBAM")
CHRSB_READS=$(samtools view -c "$OUTBAM" chrSB 2>/dev/null || echo "0")
echo "Total alignments: $TOTAL_READS"
echo "chrSB alignments: $CHRSB_READS"

echo ""
echo "Output: $OUTBAM"
echo "Date: $(date)"
echo "=== DONE ==="
