#!/bin/bash
#SBATCH --job-name=fusion_finder
#SBATCH --output=logs/fusion_finder_%j.log
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --partition=20

#
# slurm_fusion_finder.sh
#
# Run FUSION_FINDER to identify transposon-genome fusion sites
# from aligned BAM files.
#
# Prerequisites:
#   - Aligned BAMs in results/bams/
#   - input.txt configured with BAM paths
#
# Usage: sbatch scripts/slurm_fusion_finder.sh (from project root)

set -e

# === Load modules ===
module load samtools/1.17

echo "=== FUSION_FINDER Job Started ==="
echo "Date: $(date)"
echo "Node: $(hostname)"
echo "Working dir: $(pwd)"
echo ""

# Check prerequisites
if [ ! -f "input.txt" ]; then
    echo "ERROR: input.txt not found. Run sbatch from project root."
    exit 1
fi

if [ ! -f "loc_SB.txt" ]; then
    echo "ERROR: loc_SB.txt not found. Run sbatch from project root."
    exit 1
fi

# Create output directories
mkdir -p results/fusions/working

# Run fusion finder
echo "Running fusion finder on $(wc -l < input.txt) samples..."
perl scripts/fusion_finder.pl

echo ""
echo "=== DONE ==="
echo "Date: $(date)"
echo "Results: results/fusions/fusions_all.txt"
