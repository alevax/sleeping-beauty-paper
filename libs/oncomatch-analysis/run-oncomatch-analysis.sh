#!/bin/bash

# Get the directory of the script
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Change the working directory to the script directory
cd "${DIR}/../../"

mkdir -p experiments/oncomatch-analysis/

mkdir -p experiments/oncomatch-analysis/processed_data/
mkdir -p experiments/oncomatch-analysis/processed_data/gexpr/
mkdir -p experiments/oncomatch-analysis/processed_data/metadata/
mkdir -p experiments/oncomatch-analysis/processed_data/oncomatch/
mkdir -p experiments/oncomatch-analysis/processed_data/viper/

mkdir -p experiments/oncomatch-analysis/reports/

# STEP00_preprocess_matrices
echo "@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@"
echo "-------------------------------------------------------------------"
echo "----------- STEP00_SleepingBeauty_preprocess_matrices.R -----------"
echo "-------------------------------------------------------------------"
echo "@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@"
Rscript libs/oncomatch-analysis/R/STEP00_preprocess_matrices/STEP00_SleepingBeauty_preprocess_matrices.R
echo ""

echo "@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@"
echo "-------------------------------------------------------------------"
echo "------------ STEP00_SU2CEastCoast_preprocess_matrices.R -----------"
echo "-------------------------------------------------------------------"
echo "@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@"
Rscript libs/oncomatch-analysis/R/STEP00_preprocess_matrices/STEP00_SU2CEastCoast_preprocess_matrices.R
echo ""

# STEP01_preprocess_metadata
echo "@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@"
echo "-------------------------------------------------------------------"
echo "--------------- STEP01_Beltran_preprocess_metadata.R --------------"
echo "-------------------------------------------------------------------"
echo "@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@"
Rscript libs/oncomatch-analysis/R/STEP01_preprocess_metadata/STEP01_Beltran_preprocess_metadata.R
echo ""

echo "@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@"
echo "-------------------------------------------------------------------"
echo "----------- STEP01_SleepingBeauty_preprocess_metadata.R -----------"
echo "-------------------------------------------------------------------"
echo "@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@"
Rscript libs/oncomatch-analysis/R/STEP01_preprocess_metadata/STEP01_SleepingBeauty_preprocess_metadata.R
echo ""

echo "@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@"
echo "-------------------------------------------------------------------"
echo "------------ STEP01_SU2CEastCoast_preprocess_metadata.R -----------"
echo "-------------------------------------------------------------------"
echo "@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@"
Rscript libs/oncomatch-analysis/R/STEP01_preprocess_metadata/STEP01_SU2CEastCoast_preprocess_metadata.R
echo ""

echo "@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@"
echo "-------------------------------------------------------------------"
echo "------------ STEP01_SU2CWestCoast_preprocess_metadata.R -----------"
echo "-------------------------------------------------------------------"
echo "@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@"
Rscript libs/oncomatch-analysis/R/STEP01_preprocess_metadata/STEP01_SU2CWestCoast_preprocess_metadata.R
echo ""

echo "@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@"
echo "-------------------------------------------------------------------"
echo "--------------- STEP02_run_viper_using_zscore_sig.R ---------------"
echo "-------------------------------------------------------------------"
echo "@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@"
Rscript libs/oncomatch-analysis/R/STEP02_run_viper_using_zscore_sig.R
echo ""

echo "@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@"
echo "-------------------------------------------------------------------"
echo "---------------------- STEP03_run_OncoMatch.R ---------------------"
echo "-------------------------------------------------------------------"
echo "@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@"
Rscript libs/oncomatch-analysis/R/STEP03_run_OncoMatch.R
echo ""

echo "@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@"
echo "-------------------------------------------------------------------"
echo "---------------- STEP04_plot_gemm_SU2C_Beltran_om.R ---------------"
echo "-------------------------------------------------------------------"
echo "@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@"
Rscript libs/oncomatch-analysis/R/STEP04_plot_gemm_SU2C_Beltran_om.R
echo ""

echo "@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@"
echo "-------------------------------------------------------------------"
echo "------------- STEP05_plot_om_individual_and_heatmaps.R ------------"
echo "-------------------------------------------------------------------"
echo "@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@"
Rscript libs/oncomatch-analysis/R/STEP05_plot_om_individual_and_heatmaps.R
echo ""

echo "@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@"
echo "-------------------------------------------------------------------"
echo "------------------- STEP06_oncomatch_bar_plot.R -------------------"
echo "-------------------------------------------------------------------"
echo "@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@"
Rscript libs/oncomatch-analysis/R/STEP06_oncomatch_bar_plot.R
echo ""

echo "@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@"
echo "-------------------------------------------------------------------"
echo "--------------------- STEP07_oncomatch_stats.R --------------------"
echo "-------------------------------------------------------------------"
echo "@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@"
Rscript libs/oncomatch-analysis/R/STEP07_oncomatch_stats.R
echo ""
