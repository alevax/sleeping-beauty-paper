# Libraries
suppressPackageStartupMessages(library(tidyverse))

# Helper functions
getNamesReformatted <- function(name_vector, name_type){
  name_vector_unformatted <- name_vector
  name_vector_formatted <- name_vector %>%
    stringr::str_replace_all("\\.", "-") %>%
    stringr::str_replace_all("_", "-")
  my_name_df <- as_tibble(data.frame(name_vector_unformatted, name_vector_formatted))
  colnames(my_name_df) <- c(name_type, "formatted")
  return(my_name_df)
}
mergeMatrices <- function(mat1, mat2){
  shared_rownames <- intersect(rownames(mat1), rownames(mat2))
  mat1 <- mat1[shared_rownames,]
  mat2 <- mat2[shared_rownames,]
  matM <- cbind(mat1, mat2)
  return(matM)
}
getEmptyTPMMatrix <- function(SU2CEastCoast_gExpr_fpkm){
  SU2CEastCoast_gExpr_tpm <- matrix(NA,
                                    nrow = nrow(SU2CEastCoast_gExpr_fpkm),
                                    ncol = ncol(SU2CEastCoast_gExpr_fpkm))
  colnames(SU2CEastCoast_gExpr_tpm) <- colnames(SU2CEastCoast_gExpr_fpkm)
  rownames(SU2CEastCoast_gExpr_tpm) <- rownames(SU2CEastCoast_gExpr_fpkm)
  return(SU2CEastCoast_gExpr_tpm)
}
fpkmToTpm <- function(fpkm) {
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}


# Load PolyA data
SU2CEastCoast_gExpr_fpkmPolyA <- read.csv("data/SU2C_EastCoast/prad_su2c_2019/data_mrna_seq_fpkm_polya.txt", sep="\t") %>%
  as_tibble() %>%
  group_by(Hugo_Symbol) %>%
  summarise_all(.funs = sum) %>%
  ungroup() %>%
  data.frame() %>%
  column_to_rownames("Hugo_Symbol") %>%
  as.matrix()
data_clinical_sample <- as_tibble(read.csv("data/SU2C_EastCoast/prad_su2c_2019/data_clinical_sample.txt", sep = "\t"))

metadata_names_df <- getNamesReformatted(data_clinical_sample$X.Sample.Identifier, "metadata_names")
polyA_names_df <- getNamesReformatted(colnames(SU2CEastCoast_gExpr_fpkmPolyA), "fpkmPolyA_colnames")

# Filter names that are shared with the metadata
polyA_names_df <- polyA_names_df  %>%
  dplyr::filter(formatted %in% metadata_names_df$formatted)

# Switch formatted names in our names DF to match those in the original metadata names
for(i in 1:nrow(polyA_names_df)){
  metadata_name_i <- metadata_names_df %>%
    filter(formatted == polyA_names_df$formatted[i]) %>%
    pull(metadata_names)
  polyA_names_df$formatted[i] <- metadata_name_i
}
# table(polyA_names_df$formatted %in% metadata_names_df$metadata_names)
# TRUE 
# 120 
colnames(polyA_names_df)[2] <- "metadata_names"

# Select the samples that appear with the metadata
SU2CEastCoast_gExpr_fpkmPolyA <- SU2CEastCoast_gExpr_fpkmPolyA[, polyA_names_df$fpkmPolyA_colnames]
# dim(SU2CEastCoast_gExpr_fpkmPolyA)
# [1] 19270   120

# Replace the colnames to match those in the original metadata names
for(i in 1:ncol(SU2CEastCoast_gExpr_fpkmPolyA)){
  metadata_name_i <- polyA_names_df %>%
    filter(fpkmPolyA_colnames == colnames(SU2CEastCoast_gExpr_fpkmPolyA)[i]) %>%
    pull(metadata_names)
  colnames(SU2CEastCoast_gExpr_fpkmPolyA)[i] <- metadata_name_i
}

SU2CEastCoast_gExpr_tpmPolyA <- getEmptyTPMMatrix(SU2CEastCoast_gExpr_fpkmPolyA)
for(i in 1:ncol(SU2CEastCoast_gExpr_tpmPolyA)){
  SU2CEastCoast_gExpr_tpmPolyA[,i] <- fpkmToTpm(SU2CEastCoast_gExpr_fpkmPolyA[,i])
}
saveRDS(SU2CEastCoast_gExpr_tpmPolyA, "experiments/oncomatch-analysis/processed_data/gexpr/SU2CEastCoast_gExpr_polyA_tpm.rds")
