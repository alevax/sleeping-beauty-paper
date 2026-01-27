#!/bin/bash
#SBATCH --job-name=salmon_index
#SBATCH --output=logs/salmon_index_%j.log
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=20

#
# slurm_build_salmon_index.sh
#
# Build Salmon index from GENCODE mm10 transcriptome
#
# Usage: sbatch scripts/slurm_build_salmon_index.sh (from isoform_switching/)
#

set -e

# Salmon is available at /usr/local/bin/salmon (v1.8.0)

# === Configuration ===
GENCODE_VERSION="vM25"
TRANSCRIPTOME_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.transcripts.fa.gz"
GTF_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz"

REF_DIR="reference"
TRANSCRIPTOME="$REF_DIR/gencode.${GENCODE_VERSION}.transcripts.fa"
GTF="$REF_DIR/gencode.${GENCODE_VERSION}.annotation.gtf"
INDEX_DIR="$REF_DIR/salmon_index"

THREADS=${SLURM_CPUS_PER_TASK:-8}

# === Setup ===
mkdir -p "$REF_DIR"
mkdir -p logs

echo "=== Building Salmon Index ==="
echo "Date: $(date)"
echo "Node: $(hostname)"
echo "GENCODE version: $GENCODE_VERSION"
echo "Threads: $THREADS"
echo ""

# === Download transcriptome if needed ===
if [ ! -f "$TRANSCRIPTOME" ]; then
    echo "Downloading GENCODE transcriptome..."
    wget -O "${TRANSCRIPTOME}.gz" "$TRANSCRIPTOME_URL"
    gunzip "${TRANSCRIPTOME}.gz"
else
    echo "Transcriptome already exists: $TRANSCRIPTOME"
fi

# === Download GTF if needed ===
if [ ! -f "$GTF" ]; then
    echo "Downloading GENCODE GTF..."
    wget -O "${GTF}.gz" "$GTF_URL"
    gunzip "${GTF}.gz"
else
    echo "GTF already exists: $GTF"
fi

# === Build Salmon index ===
if [ -d "$INDEX_DIR" ]; then
    echo "Salmon index already exists: $INDEX_DIR"
    echo "Delete it first if you want to rebuild."
else
    echo ""
    echo "Building Salmon index..."
    salmon index \
        -t "$TRANSCRIPTOME" \
        -i "$INDEX_DIR" \
        -p $THREADS \
        --gencode

    echo ""
    echo "Salmon index built: $INDEX_DIR"
fi

echo ""
echo "Date: $(date)"
echo "=== Done ==="
