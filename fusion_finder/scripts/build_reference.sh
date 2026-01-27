#!/bin/bash
#
# build_reference.sh
#
# Build mm10+chrSB reference for fusion detection analysis
# Based on Temiz et al. (Genome Research) approach:
#   - mm10 genome with En2 gene masked (used in T2/Onc transposon)
#   - chrSB (2163bp T2/Onc2 transposon) appended as extra chromosome
#
# Prerequisites:
#   - Download mm10.fa.gz from UCSC:
#     wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz
#   - Download mm10 GTF from GENCODE:
#     wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz
#
# Required tools: samtools, bedtools, hisat2
#
# Usage: bash scripts/build_reference.sh (from project root)

set -e  # Exit on error

# === CONFIGURATION (relative paths) ===
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

REF_DIR="$PROJECT_DIR/reference"
cd "$REF_DIR"

# Input files
MM10_FA="mm10.fa"
MM10_GTF="gencode.vM25.annotation.gtf"
CHRSB_FA="chrSB.fa"

# Output files
MM10_MASKED="mm10_En2_masked.fa"
COMBINED_REF="mm10_chrSB.fa"
HISAT2_INDEX="mm10_chrSB"

# Number of threads
THREADS=4

# === STEP 1: Check prerequisites ===
echo "=== Step 1: Checking prerequisites ==="

if [ ! -f "$MM10_FA" ]; then
    echo "ERROR: $MM10_FA not found!"
    echo "Download from: https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz"
    echo "Then: gunzip mm10.fa.gz"
    exit 1
fi

if [ ! -f "$CHRSB_FA" ]; then
    echo "ERROR: $CHRSB_FA not found!"
    exit 1
fi

# Check tools
for tool in samtools bedtools hisat2-build; do
    if ! command -v $tool &> /dev/null; then
        echo "ERROR: $tool not found in PATH"
        exit 1
    fi
done

echo "All prerequisites met!"

# === STEP 2: Create En2 mask BED file ===
echo ""
echo "=== Step 2: Creating En2 mask region ==="

# En2 gene coordinates in mm10 (GRCm38)
# En2 is on chr5, the splice acceptor (SA) from En2 was used in T2/Onc
# Full En2 gene region to mask: chr5:28,148,329-28,157,139 (- strand)
# Source: UCSC Genome Browser mm10

cat > En2_mask.bed << 'EOF'
chr5	28148329	28157139	En2	0	-
EOF

echo "Created En2_mask.bed:"
cat En2_mask.bed

# === STEP 3: Mask En2 in mm10 ===
echo ""
echo "=== Step 3: Masking En2 gene in mm10 ==="

# Index mm10 if not already done
if [ ! -f "${MM10_FA}.fai" ]; then
    echo "Indexing mm10.fa..."
    samtools faidx "$MM10_FA"
fi

# Mask En2 region with Ns
echo "Running bedtools maskfasta..."
bedtools maskfasta -fi "$MM10_FA" -bed En2_mask.bed -fo "$MM10_MASKED"

echo "Created $MM10_MASKED"

# === STEP 4: Combine mm10_masked + chrSB ===
echo ""
echo "=== Step 4: Creating combined reference ==="

cat "$MM10_MASKED" "$CHRSB_FA" > "$COMBINED_REF"

echo "Created $COMBINED_REF"

# Index combined reference
echo "Indexing combined reference..."
samtools faidx "$COMBINED_REF"

# Verify chrSB is present
echo ""
echo "Verifying chrSB in combined reference:"
grep "^chrSB" "${COMBINED_REF}.fai"

# === STEP 5: Build HISAT2 index ===
echo ""
echo "=== Step 5: Building HISAT2 index ==="
echo "This will take ~30-60 minutes..."

mkdir -p hisat2_index

# Build with splice sites from GTF if available
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
    echo "Download from: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz"
    hisat2-build -p $THREADS "$COMBINED_REF" "hisat2_index/${HISAT2_INDEX}"
fi

echo ""
echo "=== DONE ==="
echo ""
echo "Reference files created:"
ls -lh "$COMBINED_REF"*
echo ""
echo "HISAT2 index files:"
ls -lh hisat2_index/

echo ""
echo "Next step: Run alignment with align_samples.sh"
