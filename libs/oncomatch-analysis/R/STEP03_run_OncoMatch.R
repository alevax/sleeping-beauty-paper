# Load libraries
suppressMessages(suppressWarnings(library(optparse)))
suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(crayon)))
suppressMessages(suppressWarnings(library(viper)))

# Functions
source("libs/tools/oncomatch_funcs.R")

# Analysis
GEMM_vpmat=readRDS("experiments/oncomatch-analysis/processed_data/viper/SleepingBeauty_viper_zscore_SU2CNets_SBNet.rds")
SU2CEastCoast_vpmat=readRDS("experiments/oncomatch-analysis/processed_data/viper/SU2CEastCoast_viper_zscore_SU2CNets.rds")
SU2CWestCoast_vpmat=readRDS("experiments/oncomatch-analysis/processed_data/viper/SU2CWestCoast_viper_zscore_SU2CNets.rds")
Beltran_vpmat=readRDS("experiments/oncomatch-analysis/processed_data/viper/Beltran_viper_zscore_SU2CNets.rds")

SU2CEastCoast_om_res <- OncoMatch(
  vpmat_to_test = GEMM_vpmat,
  vpmat_for_cMRs = SU2CEastCoast_vpmat,
  both_ways = FALSE
)
SU2CWestCoast_om_res <- OncoMatch(
  vpmat_to_test = GEMM_vpmat,
  vpmat_for_cMRs = SU2CWestCoast_vpmat,
  both_ways = FALSE
)
Beltran_om_res <- OncoMatch(
  vpmat_to_test = GEMM_vpmat,
  vpmat_for_cMRs = Beltran_vpmat,
  both_ways = FALSE
)

saveRDS(SU2CEastCoast_om_res, "experiments/oncomatch-analysis/processed_data/oncomatch/SleepingBeauty_SU2CEastCoast_OncoMatch.rds")
saveRDS(SU2CWestCoast_om_res, "experiments/oncomatch-analysis/processed_data/oncomatch/SleepingBeauty_SU2CWestCoast_OncoMatch.rds")
saveRDS(Beltran_om_res, "experiments/oncomatch-analysis/processed_data/oncomatch/SleepingBeauty_Beltran_OncoMatch.rds")
