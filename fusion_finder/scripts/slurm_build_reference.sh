#!/bin/bash
#SBATCH --job-name=build_mm10_chrSB
#SBATCH --output=logs/build_reference_%j.log
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=256G
#SBATCH --partition=20

#
# slurm_build_reference.sh
#
# Build mm10+chrSB reference for fusion detection analysis
# Based on Temiz et al. (Genome Research) approach
#
# Usage: sbatch scripts/slurm_build_reference.sh (from project root)
#

set -e  # Exit on error

THREADS=${SLURM_CPUS_PER_TASK:-16}

# === Load modules ===
module load samtools/1.17
module load bedtools2/2.29.2
module load hisat2/2.2.1

cd reference

echo "=== Build Reference Job Started ==="
echo "Date: $(date)"
echo "Node: $(hostname)"
echo "Working dir: $(pwd)"
echo "Threads: $THREADS"
echo ""

# Input files
MM10_FA="mm10.fa"
MM10_GTF="gencode.vM25.annotation.gtf"
CHRSB_FA="chrSB.fa"

# Output files
MM10_MASKED="mm10_En2_masked.fa"
COMBINED_REF="mm10_chrSB.fa"
HISAT2_INDEX="mm10_chrSB"

# === STEP 1: Check prerequisites ===
echo "=== Step 1: Checking prerequisites ==="

if [ ! -f "$MM10_FA" ]; then
    echo "Downloading mm10.fa from UCSC..."
    wget -q https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz
    gunzip mm10.fa.gz
fi

if [ ! -f "$MM10_GTF" ]; then
    echo "Downloading GENCODE mm10 GTF..."
    wget -q https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz
    gunzip gencode.vM25.annotation.gtf.gz
fi

if [ ! -f "$CHRSB_FA" ]; then
    echo "ERROR: $CHRSB_FA not found!"
    exit 1
fi

echo "All prerequisites met!"

# === STEP 2: Create En2 mask BED file ===
echo ""
echo "=== Step 2: Creating En2 mask region ==="

# En2 gene coordinates in mm10 (GRCm38)
# Full En2 gene region to mask: chr5:28,148,329-28,157,139 (- strand)
cat > En2_mask.bed << 'EOF'
chr5	28148329	28157139	En2	0	-
EOF

echo "Created En2_mask.bed"

# === STEP 3: Mask En2 in mm10 ===
echo ""
echo "=== Step 3: Masking En2 gene in mm10 ==="

if [ ! -f "${MM10_FA}.fai" ]; then
    echo "Indexing mm10.fa..."
    samtools faidx "$MM10_FA"
fi

echo "Running bedtools maskfasta..."
bedtools maskfasta -fi "$MM10_FA" -bed En2_mask.bed -fo "$MM10_MASKED"
echo "Created $MM10_MASKED"

# === STEP 4: Combine mm10_masked + chrSB ===
echo ""
echo "=== Step 4: Creating combined reference ==="

cat "$MM10_MASKED" "$CHRSB_FA" > "$COMBINED_REF"
echo "Created $COMBINED_REF"

samtools faidx "$COMBINED_REF"

echo "Verifying chrSB in combined reference:"
grep "^chrSB" "${COMBINED_REF}.fai"

# === STEP 5: Build HISAT2 index ===
echo ""
echo "=== Step 5: Building HISAT2 index ==="

mkdir -p hisat2_index

if [ -f "$MM10_GTF" ]; then
    echo "Extracting splice sites from GTF..."
    hisat2_extract_splice_sites.py "$MM10_GTF" > hisat2_index/splice_sites.txt
    hisat2_extract_exons.py "$MM10_GTF" > hisat2_index/exons.txt

    echo "Building HISAT2 index with splice site annotation..."
    hisat2-build -p $THREADS \
        --ss hisat2_index/splice_sites.txt \
        --exon hisat2_index/exons.txt \
        "$COMBINED_REF" \
        "hisat2_index/${HISAT2_INDEX}"
else
    echo "WARNING: GTF not found, building index without splice sites"
    hisat2-build -p $THREADS "$COMBINED_REF" "hisat2_index/${HISAT2_INDEX}"
fi

# === Done ===
echo ""
echo "=== DONE ==="
echo "Date: $(date)"
echo ""
echo "Reference files created:"
ls -lh "$COMBINED_REF"*
echo ""
echo "HISAT2 index files:"
ls -lh hisat2_index/
echo ""
echo "Next step: sbatch scripts/slurm_align_array.sh"
