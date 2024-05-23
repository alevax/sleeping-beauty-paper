# Load libraries
suppressMessages(suppressWarnings(library(optparse)))
suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(crayon)))
suppressMessages(suppressWarnings(library(viper)))

# Functions
source("libs/tools/oncomatch_funcs.R")
saveLargeMatToExcel <- function(mat, filename){
  mat_xlsx <- data.frame(cbind(rownames(mat), mat))
  colnames(mat_xlsx)[1] <- "rowname"
  writexl::write_xlsx(mat_xlsx, filename, col_names = TRUE)
}

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

saveLargeMatToExcel(
  SU2CEastCoast_om_res,
  "experiments/oncomatch-analysis/processed_data/oncomatch/SleepingBeauty_SU2CEastCoast_OncoMatch.xlsx"
)
saveLargeMatToExcel(
  Beltran_om_res,
  "experiments/oncomatch-analysis/processed_data/oncomatch/SleepingBeauty_Beltran_OncoMatch.xlsx"
)

SU2CEastCoast_om_res_gsea <- OncoMatch(
  vpmat_to_test = GEMM_vpmat,
  vpmat_for_cMRs = SU2CEastCoast_vpmat,
  both_ways = FALSE,
  enrichment = "gsea",
  n_per = 500
)
SU2CWestCoast_om_res_gsea <- OncoMatch(
  vpmat_to_test = GEMM_vpmat,
  vpmat_for_cMRs = SU2CWestCoast_vpmat,
  both_ways = FALSE,
  enrichment = "gsea",
  n_per = 500
)
Beltran_om_res_gsea <- OncoMatch(
  vpmat_to_test = GEMM_vpmat,
  vpmat_for_cMRs = Beltran_vpmat,
  both_ways = FALSE,
  enrichment = "gsea",
  n_per = 500
)

saveRDS(SU2CEastCoast_om_res_gsea,
        "experiments/oncomatch-analysis/processed_data/oncomatch/SleepingBeauty_SU2CEastCoast_OncoMatch_gsea.rds")
saveRDS(SU2CWestCoast_om_res_gsea,
        "experiments/oncomatch-analysis/processed_data/oncomatch/SleepingBeauty_SU2CWestCoast_OncoMatch_gsea.rds")
saveRDS(Beltran_om_res_gsea,
        "experiments/oncomatch-analysis/processed_data/oncomatch/SleepingBeauty_Beltran_OncoMatch_gsea.rds")

saveLargeMatToExcel(
  SU2CEastCoast_om_res_gsea,
  "experiments/oncomatch-analysis/processed_data/oncomatch/SleepingBeauty_SU2CEastCoast_OncoMatch_gsea.xlsx"
)
saveLargeMatToExcel(
  Beltran_om_res_gsea,
  "experiments/oncomatch-analysis/processed_data/oncomatch/SleepingBeauty_Beltran_OncoMatch_gsea.xlsx"
)




vpmat_to_test <- GEMM_vpmat
vpmat_for_cMRs <- Beltran_vpmat
tcm_size <- 50

cMRs_regulon <- generateRegulonObjectFromProteinActivityMatrix(vpmat_for_cMRs,
                                                             n_top = tcm_size/2,
                                                             justWithOnes = T)
set.seed(0)
pdf("experiments/oncomatch-analysis/reports/aREA/CMZ103_WCMC0_6_N_aREA_51.8269.pdf")
gsea( vpmat_to_test[,"CMZ103"], cMRs_regulon$WCMC0_6_N$tfmode, twoTails=TRUE, per = 0, pout = TRUE)
dev.off()

set.seed(0)
pdf("experiments/oncomatch-analysis/reports/aREA/CMZ181_WCMC4240_1_N_aREA_0.pdf")
gsea( vpmat_to_test[,"CMZ181"], cMRs_regulon$WCMC4240_1_N$tfmode, twoTails=TRUE, per = 0, pout = TRUE)
dev.off()




