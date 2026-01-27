#!/bin/bash
#
# align_samples.sh
#
# Align RNA-seq FASTQs to mm10+chrSB reference using HISAT2
# Handles both single-end and paired-end samples
#
# Prerequisites:
#   - Run build_reference.sh first
#   - HISAT2 and samtools in PATH
#
# Usage: bash scripts/align_samples.sh (from project root)

set -e  # Exit on error

# === CONFIGURATION (relative paths) ===
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

FASTQ_DIR="$PROJECT_DIR/data/fastq"
OUTPUT_DIR="$PROJECT_DIR/results/bams"
STATS_DIR="$OUTPUT_DIR/stats"
HISAT2_INDEX="$PROJECT_DIR/reference/hisat2_index/mm10_chrSB"
SAMPLES_FILE="$PROJECT_DIR/data/fastq/samples.txt"

# Number of threads
THREADS=4

# === Setup ===
mkdir -p "$OUTPUT_DIR"
mkdir -p "$STATS_DIR"

# Check HISAT2 index exists
if [ ! -f "${HISAT2_INDEX}.1.ht2" ]; then
    echo "ERROR: HISAT2 index not found at $HISAT2_INDEX"
    echo "Run build_reference.sh first!"
    exit 1
fi

echo "=== RNA-seq Alignment Pipeline ==="
echo "FASTQ dir: $FASTQ_DIR"
echo "Output dir: $OUTPUT_DIR"
echo "HISAT2 index: $HISAT2_INDEX"
echo ""

# === Alignment function for single-end ===
align_single_end() {
    local sample=$1
    local fastq=$(ls "$FASTQ_DIR"/${sample}_*R1*.fastq.gz 2>/dev/null | head -1)

    if [ -z "$fastq" ]; then
        echo "WARNING: No FASTQ found for $sample, skipping..."
        return
    fi

    local outbam="$OUTPUT_DIR/${sample}_sorted.bam"

    if [ -f "$outbam" ]; then
        echo "SKIP: $outbam already exists"
        return
    fi

    echo "Aligning $sample (single-end)..."
    echo "  Input: $fastq"

    hisat2 -p $THREADS \
        -x "$HISAT2_INDEX" \
        -U "$fastq" \
        --dta \
        --summary-file "$STATS_DIR/${sample}_alignment_summary.txt" \
        2> "$STATS_DIR/${sample}_hisat2.log" \
    | samtools sort -@ $THREADS -o "$outbam" -

    samtools index "$outbam"

    echo "  Output: $outbam"
    echo "  $(samtools view -c $outbam) total alignments"
}

# === Alignment function for paired-end ===
align_paired_end() {
    local sample=$1
    local fastq_r1=$(ls "$FASTQ_DIR"/${sample}_*R1*.fastq.gz 2>/dev/null | head -1)
    local fastq_r2=$(ls "$FASTQ_DIR"/${sample}_*R2*.fastq.gz 2>/dev/null | head -1)

    if [ -z "$fastq_r1" ] || [ -z "$fastq_r2" ]; then
        echo "WARNING: Missing FASTQ pair for $sample, skipping..."
        return
    fi

    local outbam="$OUTPUT_DIR/${sample}_sorted.bam"

    if [ -f "$outbam" ]; then
        echo "SKIP: $outbam already exists"
        return
    fi

    echo "Aligning $sample (paired-end)..."
    echo "  R1: $fastq_r1"
    echo "  R2: $fastq_r2"

    hisat2 -p $THREADS \
        -x "$HISAT2_INDEX" \
        -1 "$fastq_r1" \
        -2 "$fastq_r2" \
        --dta \
        --summary-file "$STATS_DIR/${sample}_alignment_summary.txt" \
        2> "$STATS_DIR/${sample}_hisat2.log" \
    | samtools sort -@ $THREADS -o "$outbam" -

    samtools index "$outbam"

    echo "  Output: $outbam"
    echo "  $(samtools view -c $outbam) total alignments"
}

# === Run alignments from samples.txt ===
while IFS=',' read -r sample seq_type; do
    if [ "$seq_type" == "single" ]; then
        align_single_end "$sample"
    else
        align_paired_end "$sample"
    fi
done < "$SAMPLES_FILE"

# === Summary ===
echo ""
echo "=== Alignment Summary ==="
echo ""
echo "Sample alignments completed:"
ls -lh "$OUTPUT_DIR"/*.bam 2>/dev/null | wc -l
echo ""

# Check for chrSB reads in each BAM
echo "chrSB read counts per sample:"
for bam in "$OUTPUT_DIR"/*.bam; do
    sample=$(basename "$bam" _sorted.bam)
    chrSB_reads=$(samtools view -c "$bam" chrSB 2>/dev/null || echo "0")
    echo "  $sample: $chrSB_reads reads on chrSB"
done

echo ""
echo "=== DONE ==="
echo "BAM files with chrSB alignments are in: $OUTPUT_DIR"
