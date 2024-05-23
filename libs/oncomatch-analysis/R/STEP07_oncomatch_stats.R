library(tidyverse)
getBeltranNEAnnotation <- function(){
  Beltran_metadata <- readRDS("experiments/oncomatch-analysis/processed_data/metadata/Beltran_metadata_forOMPlot.rds")
  Beltran_metadata_NE_feats <- Beltran_metadata$Neuroendocrine_Features
  names(Beltran_metadata_NE_feats) <- Beltran_metadata$Sample_Identifier
  return(Beltran_metadata_NE_feats)
}
getSU2CEastCoastNEAnnotation <- function(){
  SU2CEastCoast_metadata <- readRDS("experiments/oncomatch-analysis/processed_data/metadata/SU2C_EastCoast_metadata_forOMPlot.rds")
  SU2CEastCoast_metadata <- SU2CEastCoast_metadata %>%
    filter(!is.na(Neuroendocrine_Features)) #%>%
  SU2CEastCoast_metadata_NE_feats <- SU2CEastCoast_metadata$Neuroendocrine_Features
  names(SU2CEastCoast_metadata_NE_feats) <- SU2CEastCoast_metadata$Sample_Identifier
  return(SU2CEastCoast_metadata_NE_feats)
}
getPercHumanToGEMMMatches <- function(Human_om, SB_metadata_NE_feats, Human_metadata_NE_feats, om_threshold){
  round(sum(apply(Human_om[SB_metadata_NE_feats == "NEPC", Human_metadata_NE_feats == "Yes"] >= om_threshold, 2, any))/sum(Human_metadata_NE_feats == "Yes")*100,1)
}
getPercGEMMToHumanMatches <- function(Human_om, SB_metadata_NE_feats, Human_metadata_NE_feats, om_threshold){
  round(sum(apply(Human_om[SB_metadata_NE_feats == "NEPC", Human_metadata_NE_feats == "Yes"] >= om_threshold, 1, any))/sum(SB_metadata_NE_feats == "NEPC")*100,1)
}
run_analysis <- function(Beltran_om, Su2c_EC_om, gemm_metadata, reports_results_dir, om_threshold = 10, enrichment_method){
  Beltran_metadata_NE_feats <- getBeltranNEAnnotation()
  Su2c_EC_metadata_NE_feats <- getSU2CEastCoastNEAnnotation()
  SB_metadata_NE_feats <- gemm_metadata$Phenotype
  
  Beltran_NE_SB_NE_perc_matches <- getPercHumanToGEMMMatches(
    Beltran_om,
    SB_metadata_NE_feats,
    Beltran_metadata_NE_feats,
    om_threshold
  )
  Su2c_EC_NE_SB_NE_perc_matches <- getPercHumanToGEMMMatches(
    Su2c_EC_om,
    SB_metadata_NE_feats,
    Su2c_EC_metadata_NE_feats,
    om_threshold
  )
  
  SB_NE_Beltran_NE_perc_matches <- getPercGEMMToHumanMatches(
    Beltran_om,
    SB_metadata_NE_feats,
    Beltran_metadata_NE_feats,
    om_threshold
  )
  SB_NE_Su2c_EC_NE_perc_matches <- getPercGEMMToHumanMatches(
    Su2c_EC_om,
    SB_metadata_NE_feats,
    Su2c_EC_metadata_NE_feats,
    om_threshold
  )
  
  res_df <- data.frame(t(matrix(
    c("Beltran", "SB", Beltran_NE_SB_NE_perc_matches,
      "SU2C_EC", "SB", Su2c_EC_NE_SB_NE_perc_matches,
      "SB", "Beltran", SB_NE_Beltran_NE_perc_matches,
      "SB", "SU2C_EC", SB_NE_Su2c_EC_NE_perc_matches),
    nrow = 3
  )))
  
  colnames(res_df) <- c("Group", "Matching_To", paste0("Perc_matched (om_threshold=", om_threshold, ")"))
  writexl::write_xlsx(res_df, paste0(reports_results_dir, enrichment_method, "_match_rates.xlsx"), col_names = TRUE)
}

Beltran_om <-
  readRDS("experiments/oncomatch-analysis/processed_data/oncomatch/SleepingBeauty_Beltran_OncoMatch.rds")
Su2c_EC_om <-
  readRDS("experiments/oncomatch-analysis/processed_data/oncomatch/SleepingBeauty_SU2CEastCoast_OncoMatch.rds")
gemm_metadata <-
  readRDS("experiments/oncomatch-analysis/processed_data/metadata/SleepingBeauty_metadata_forOMPlot.rds")

Beltran_om_gsea <-
  readRDS("experiments/oncomatch-analysis/processed_data/oncomatch/SleepingBeauty_Beltran_OncoMatch_gsea.rds")
Su2c_EC_om_gsea <-
  readRDS("experiments/oncomatch-analysis/processed_data/oncomatch/SleepingBeauty_SU2CEastCoast_OncoMatch_gsea.rds")


reports_results_dir_aREA <- "experiments/oncomatch-analysis/reports/aREA/"
dir.create(reports_results_dir_aREA)
reports_results_dir_GSEA <- "experiments/oncomatch-analysis/reports/GSEA/"
dir.create(reports_results_dir_GSEA)

run_analysis(
  Beltran_om,
  Su2c_EC_om,
  gemm_metadata,
  reports_results_dir_aREA,
  om_threshold = 10,
  enrichment_method = 'aREA'
)
run_analysis(
  Beltran_om_gsea,
  Su2c_EC_om_gsea,
  gemm_metadata,
  reports_results_dir_GSEA,
  om_threshold = 5,
  enrichment_method = 'GSEA'
)

