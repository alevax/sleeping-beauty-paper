# Libraries
suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(viper)))

# Functions
GESTransform <- function(dat.mat, na.rm = TRUE){
  # generate GES
  ges.mat <- t(apply(dat.mat, 1, function(x) {
    (x - mean(x, na.rm = na.rm)) / sd(x, na.rm = na.rm)
  }))
  ges.mat
}
filterOutNPSamples <- function(SleepingBeauty_mat, SleepingBeauty_metadata){
  NPp53_sample_names <- SleepingBeauty_metadata %>%
    filter(genotype == "NPp53") %>%
    pull(Sample_Identifier)
  SleepingBeauty_mat <- SleepingBeauty_mat[, NPp53_sample_names]
  return(SleepingBeauty_mat)
}

# Load Matrices
SU2CWestCoast_gExpr_tpm <- readRDS("data/SU2C_WestCoast/SU2CWestCoast_merged_tpm_geneNames.rds")
SU2CEastCoast_gExpr_tpm <- readRDS("experiments/oncomatch-analysis/processed_data/gexpr/SU2CEastCoast_gExpr_polyA_tpm.rds")
Beltran_gExpr_tpm <- readRDS("data/Beltran-2016/beltran-sym-tpm-original.rds")
SleepingBeauty_gExpr_tpm <- readRDS("experiments/oncomatch-analysis/processed_data/gexpr/SleepingBeauty_gExpr_tpm_humanNames.rds")
# Load Networks
su2cEC_regulon <- readRDS("data/SU2C_EastCoast/SU2CEastCoast_net_pruned.rds")
su2cWC_regulon <- readRDS("data/SU2C_WestCoast/SU2CWestCoast_net_pruned.rds")
SleepingBeauty_regulon <- readRDS("data/sb-networks/SleepingBeauty_net_pruned_humanNames.rds")
# Load Metadata
SleepingBeauty_metadata <- readRDS("experiments/oncomatch-analysis/processed_data/metadata/SleepingBeauty_metadata_forOMPlot.rds")

# Analysis
# Scale transforming matrices
SleepingBeauty_gExpr_zscore <- GESTransform(SleepingBeauty_gExpr_tpm)
SU2CEastCoast_gExpr_zscore <- GESTransform(SU2CEastCoast_gExpr_tpm)
SU2CWestCoast_gExpr_zscore <- GESTransform(SU2CWestCoast_gExpr_tpm)
Beltran_gExpr_zscore <- GESTransform(Beltran_gExpr_tpm)

# Filtering out gemm NP samples
SleepingBeauty_gExpr_zscore <- filterOutNPSamples(SleepingBeauty_gExpr_zscore, SleepingBeauty_metadata)

# Compute the VIPER
SleepingBeauty_viper <- viper(SleepingBeauty_gExpr_zscore,
                              # list(SleepingBeauty_regulon, su2cEC_regulon, su2cWC_regulon),
                              list(SleepingBeauty_regulon, su2cEC_regulon),
                              method = 'none',
                              eset.filter = FALSE,
                              mvws = 10)
SU2CEastCoast_viper <- viper(SU2CEastCoast_gExpr_zscore,
                             # list(su2cEC_regulon, su2cWC_regulon),
                             su2cEC_regulon,
                             method = 'none',
                             eset.filter = FALSE,
                             mvws = 10)
SU2CWestCoast_viper <- viper(SU2CWestCoast_gExpr_zscore,
                             # list(su2cEC_regulon, su2cWC_regulon),
                             su2cEC_regulon,
                             method = 'none',
                             eset.filter = FALSE,
                             mvws = 10)
Beltran_viper <- viper(Beltran_gExpr_zscore,
                       # list(su2cEC_regulon, su2cWC_regulon),
                       su2cEC_regulon,
                       method = 'none',
                       eset.filter = FALSE,
                       mvws = 10)

saveRDS(SleepingBeauty_viper, "experiments/oncomatch-analysis/processed_data/viper/SleepingBeauty_viper_zscore_SU2CNets_SBNet.rds")
saveRDS(SU2CEastCoast_viper, "experiments/oncomatch-analysis/processed_data/viper/SU2CEastCoast_viper_zscore_SU2CNets.rds")
saveRDS(SU2CWestCoast_viper, "experiments/oncomatch-analysis/processed_data/viper/SU2CWestCoast_viper_zscore_SU2CNets.rds")
saveRDS(Beltran_viper, "experiments/oncomatch-analysis/processed_data/viper/Beltran_viper_zscore_SU2CNets.rds")
