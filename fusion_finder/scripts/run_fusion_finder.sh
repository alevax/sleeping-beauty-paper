#!/bin/bash
#
# run_fusion_finder.sh
#
# Run FUSION_FINDER to identify transposon-genome fusion sites
# from aligned BAM files.
#
# Prerequisites:
#   - Aligned BAMs in results/bams/
#   - input.txt configured with BAM paths
#   - samtools in PATH
#
# Usage: bash scripts/run_fusion_finder.sh (from project root)

set -e

# === CONFIGURATION (relative paths) ===
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

cd "$PROJECT_DIR"

# Check prerequisites
if [ ! -f "input.txt" ]; then
    echo "ERROR: input.txt not found in project root"
    exit 1
fi

if [ ! -f "loc_SB.txt" ]; then
    echo "ERROR: loc_SB.txt not found in project root"
    exit 1
fi

if ! command -v samtools &> /dev/null; then
    echo "ERROR: samtools not found in PATH"
    exit 1
fi

echo "=== FUSION_FINDER ==="
echo "Date: $(date)"
echo "Project dir: $PROJECT_DIR"
echo ""

# Create output directories
mkdir -p results/fusions/working

# Run fusion finder
echo "Running fusion finder..."
perl scripts/fusion_finder.pl

echo ""
echo "=== DONE ==="
echo "Results: results/fusions/fusions_all.txt"
