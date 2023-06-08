#!/bin/bash

# Get the directory of the script
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# Change the working directory to the script directory
cd "${DIR}"

mkdir -p experiments/

echo "@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@"
echo "-------------------------------------------------------------------"
echo "----------------- ALL SAMPLES SB RNA SEQ ANALYSIS -----------------"
echo "-------------------------------------------------------------------"
echo "@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@"
mkdir -p experiments/all-samples-sb-rna-seq-analysis/
mkdir -p experiments/all-samples-sb-rna-seq-analysis/processed_data/
mkdir -p experiments/all-samples-sb-rna-seq-analysis/reports/
# Rscript libs/all-samples-sb-rna-seq-analysis-redo.R \
Rscript libs/all-samples-sb-rna-seq-analysis.R \
> experiments/all-samples-sb-rna-seq-analysis/all-samples-sb-rna-seq-analysis-log.txt 2>&1
echo ""

echo "@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@"
echo "-------------------------------------------------------------------"
echo "----------------------- INTEGRATIVE ANALYSIS ----------------------"
echo "-------------------------------------------------------------------"
echo "@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@"
mkdir -p experiments/integrative-analysis/
mkdir -p experiments/integrative-analysis/processed_data/
mkdir -p experiments/integrative-analysis/reports/
Rscript libs/integrative_analysis.R \
> experiments/integrative-analysis/integrative-analysis-log.txt 2>&1
echo ""

echo "@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@"
echo "-------------------------------------------------------------------"
echo "----------------------- SB FIGURE GENERATOR -----------------------"
echo "-------------------------------------------------------------------"
echo "@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@"
mkdir -p experiments/sb-figure-generator/
mkdir -p experiments/sb-figure-generator/processed_data/
mkdir -p experiments/sb-figure-generator/reports/
Rscript libs/sb-figure-generator.R \
> experiments/sb-figure-generator/sb-figure-generator-log.txt 2>&1
echo ""

echo "@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@"
echo "-------------------------------------------------------------------"
echo "------------------------ ONCOMATCH ANALYSIS -----------------------"
echo "-------------------------------------------------------------------"
echo "@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@"
mkdir -p experiments/survival-analysis/
mkdir -p experiments/survival-analysis/processed_data/
mkdir -p experiments/survival-analysis/reports/

mkdir -p experiments/oncomatch-analysis/
mkdir -p experiments/oncomatch-analysis/processed_data/
mkdir -p experiments/oncomatch-analysis/reports/
bash libs/oncomatch-analysis/run-oncomatch-analysis.sh \
> experiments/oncomatch-analysis/oncomatch-analysis-log.txt 2>&1
echo ""

echo "@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@"
echo "-------------------------------------------------------------------"
echo "------------------------- NETWORK ANALYSIS ------------------------"
echo "-------------------------------------------------------------------"
echo "@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@"
mkdir -p experiments/network-analysis/
mkdir -p experiments/network-analysis/processed_data/
mkdir -p experiments/network-analysis/reports/
Rscript libs/network_analysis.R \
> experiments/network-analysis/network-analysis-log.txt 2>&1
echo ""

echo "@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@"
echo "-------------------------------------------------------------------"
echo "------------------------- NEPC SIG ANALYSIS -----------------------"
echo "-------------------------------------------------------------------"
echo "@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@"
mkdir -p experiments/NEPC-sig-analysis/
mkdir -p experiments/NEPC-sig-analysis/processed_data/
mkdir -p experiments/NEPC-sig-analysis/reports/
Rscript libs/NEPC_signature_analysis.R \
> experiments/NEPC-sig-analysis/NEPC-sig-analysis-log.txt 2>&1
echo ""

echo "@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@"
echo "-------------------------------------------------------------------"
echo "------------------------- SURVIVAL ANALYSIS -----------------------"
echo "-------------------------------------------------------------------"
echo "@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@"
Rscript libs/survival_analysis_complete.R \
> experiments/survival-analysis/survival-analysis-log.txt 2>&1
echo ""
