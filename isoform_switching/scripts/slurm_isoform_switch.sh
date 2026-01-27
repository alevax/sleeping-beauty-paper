#!/bin/bash
#SBATCH --job-name=isoform_switch
#SBATCH --output=logs/isoform_switch_%j.log
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --partition=20

#
# slurm_isoform_switch.sh
#
# Run IsoformSwitchAnalyzeR analysis
# This is compute-intensive: DEXSeq testing, ORF analysis, consequence prediction
#
# Usage: sbatch scripts/slurm_isoform_switch.sh (from isoform_switching/)
#

set -e

echo "=== Isoform Switch Analysis ==="
echo "Date: $(date)"
echo "Node: $(hostname)"
echo "Job ID: $SLURM_JOB_ID"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Memory: $SLURM_MEM_PER_NODE MB"
echo ""

# Load R modules
module load R/4.2
module load R-site-library/4.2/2025q1

echo "R version:"
R --version | head -1
echo ""

# Change to working directory
cd /lab/barcheese01/mdiberna/sleeping_beauty/isoform_switching

echo "Working directory: $(pwd)"
echo ""

# Run the analysis
echo "Starting IsoformSwitchAnalyzeR analysis..."
echo ""

Rscript R/isoform_switch_analysis.R

echo ""
echo "Date: $(date)"
echo "=== Done ==="
