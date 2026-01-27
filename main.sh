#!/bin/bash

# Get the directory of the script
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# Change the working directory to the script directory
cd "${DIR}"

mkdir -p experiments/

echo "@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@"
echo "-------------------------------------------------------------------"
echo "------------------------- GENERATE TABLE --------------------------"
echo "-------------------------------------------------------------------"
echo "@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@"
mkdir -p experiments/table/
Rscript libs/table.R \
> experiments/table/table-analysis-log.txt 2>&1
echo ""

echo "@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@"
echo "-------------------------------------------------------------------"
echo "-------------------------------------------------------------------"
echo "------------------------- CO-OCCURRENCE ANALYSIS ------------------"
echo "-------------------------------------------------------------------"
echo "@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@"
Rscript libs/SB-co-occurrence.R \
> experiments/table/co-occurrence-analysis-log.txt 2>&1
echo ""

echo "@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@"
echo "-------------------------------------------------------------------"
echo "-------------------------------------------------------------------"
echo "------------------------- CIRCOS ANALYSIS -------------------------"
echo "-------------------------------------------------------------------"
echo "@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@"
mkdir -p experiments/circos/
Rscript libs/circos.R \
> experiments/circos/circos-analysis-log.txt 2>&1
echo ""

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

echo "@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@"
echo "-------------------------------------------------------------------"
echo "------------------------- LOLLIPOP ANALYSIS -----------------------"
echo "-------------------------------------------------------------------"
echo "@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@"
mkdir -p experiments/lollipop/
Rscript libs/lollipop.R \
> experiments/lollipop/lollipop-analysis-log.txt 2>&1
echo ""

echo "@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@"
echo "-------------------------------------------------------------------"
echo "-------------------- FUSION FINDER ANALYSIS -----------------------"
echo "-------------------------------------------------------------------"
echo "@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@"
Rscript libs/fusion_finder_analysis.R \
> experiments/fusion-finder/fusion-finder-analysis-log.txt 2>&1
echo ""

echo "@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@"
echo "-------------------------------------------------------------------"
echo "------------------ ISOFORM SWITCHING ANALYSIS ---------------------"
echo "-------------------------------------------------------------------"
echo "@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@"
mkdir -p experiments/isoform-switching/
mkdir -p experiments/isoform-switching/processed_data/
mkdir -p experiments/isoform-switching/reports/
Rscript libs/isoform-switching/01_cis_isoform_switch_overlap.R \
> experiments/isoform-switching/01-cis-overlap-log.txt 2>&1
Rscript libs/isoform-switching/02_cis_switch_consequences.R \
> experiments/isoform-switching/02-consequences-log.txt 2>&1
echo ""

echo "@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@"
echo "-------------------------------------------------------------------"
echo "------------------ CHIMERIC ISOFORMS ANALYSIS ---------------------"
echo "-------------------------------------------------------------------"
echo "@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@"
Rscript libs/chimeric_isoforms_analysis.R \
> experiments/chimeric-isoforms/chimeric-isoforms-analysis-log.txt 2>&1
echo ""
