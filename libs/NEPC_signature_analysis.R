# Libraries
library(tidyverse)
library(limma)

# Functions
getTopCISGenes <- function(){
  top_hits_tbl <- readxl::read_xlsx("experiments/integrative-analysis/processed_data/top-hits-for-validation-final-fisher-integration-table.xlsx")
  cis_genes_top <- top_hits_tbl$gene_human
  return(cis_genes_top)
}
getAllCISGenes <- function(){
  complete_tbl <- readxl::read_xlsx("experiments/integrative-analysis/reports/nepc-signatures-table-final-fisher-integration-table.xlsx")
  cis_genes_all <- complete_tbl %>%
    filter(is_cis_gene=="Yes") %>%
    pull(gene_human) %>%
    na.omit() %>%
    as.character()
  return(cis_genes_all)
}
# Form an appropriate design matrix for the two RNA sources and fit linear models. The design matrix
# has two columns. The first represents log-expression in the wild-type and the second represents the
# log-ratio between the mutant and wild-type samples. See Section 9.2 for more details on the design
# matrix.
# > design <- cbind(WT=1, MUvsWT=targets$Genotype=="mu")
runLimma <- function(gExpr_tpm, NE_feats){
  stopifnot(identical(colnames(gExpr_tpm), names(NE_feats)))
  design <- data.frame("NE" = as.integer(NE_feats == "Yes"),
                       "Adeno" = as.integer(NE_feats == "No"))
  fit <- limma::lmFit(gExpr_tpm <- log2(gExpr_tpm + 1), design)
  contrast.matrix <- makeContrasts(NEvsAdeno=NE-Adeno, levels=design)
  fit2 <- limma::contrasts.fit(fit, contrast.matrix)
  fit2 <- limma::eBayes(fit2, trend=TRUE)
  res <- limma::topTable(fit2, number = nrow(gExpr_tpm))
  res <- res %>%
    rownames_to_column() %>%
    as_tibble()
  return(res)
}
fishers_method <- function(p_values) {
  neg_log_p <- -log(p_values)
  fisher_statistic <- -2 * sum(neg_log_p)
  # combined_p <- pchisq(fisher_statistic, df = 2 * length(p_values), lower.tail = TRUE)
  log_combined_p <- pchisq(fisher_statistic, df = 2 * length(p_values), lower.tail = TRUE, log.p = TRUE)
  
  combined_p <- as.numeric(sprintf("%.2e", exp(log_combined_p)))
  return(combined_p)
}
fishers_method_for_mutate <- function(p1,p2,p3) {
  p_values <- c(p1,p2,p3)
  chi_squared_statistic <- -2 * sum(log(p_values))
  combined_p <- 1 - pchisq(chi_squared_statistic, df = 2 * length(p_values))
  return(combined_p)
}
getIntegratedGExprSigDFCISGenes <- function(Beltran_NEvsAdeno_gExpr_sig_df_cisGenes,
                                            SU2CEastCoast_NEvsAdeno_gExpr_sig_df_cisGenes,
                                            SU2CWestCoast_NEvsAdeno_gExpr_sig_df_cisGenes,
                                            combined_gExpr_sig_df_cisGenes_All){
  Beltran_NEvsAdeno_gExpr_sig_df_cisGenes$rowname <- NULL
  Beltran_NEvsAdeno_gExpr_sig_df_cisGenes$Dataset <- "Beltran"
  SU2CEastCoast_NEvsAdeno_gExpr_sig_df_cisGenes$Dataset <- "SU2CEastCoast"
  SU2CWestCoast_NEvsAdeno_gExpr_sig_df_cisGenes$Dataset <- "SU2CWestCoast"
  integratedPValueDF <- dplyr::select(combined_gExpr_sig_df_cisGenes_All, ID, fischer_integrated_adj.P.Val)
  Beltran_NEvsAdeno_gExpr_sig_df_cisGenes <- dplyr::left_join(
    Beltran_NEvsAdeno_gExpr_sig_df_cisGenes,
    integratedPValueDF,
    by = c("ID" = "ID")
  )
  SU2CEastCoast_NEvsAdeno_gExpr_sig_df_cisGenes <- dplyr::left_join(
    SU2CEastCoast_NEvsAdeno_gExpr_sig_df_cisGenes,
    integratedPValueDF,
    by = c("ID" = "ID")
  )
  SU2CWestCoast_NEvsAdeno_gExpr_sig_df_cisGenes <- dplyr::left_join(
    SU2CWestCoast_NEvsAdeno_gExpr_sig_df_cisGenes,
    integratedPValueDF,
    by = c("ID" = "ID")
  )
  stopifnot(identical(colnames(Beltran_NEvsAdeno_gExpr_sig_df_cisGenes),
                      colnames(SU2CEastCoast_NEvsAdeno_gExpr_sig_df_cisGenes)))
  stopifnot(identical(colnames(Beltran_NEvsAdeno_gExpr_sig_df_cisGenes),
                      colnames(SU2CWestCoast_NEvsAdeno_gExpr_sig_df_cisGenes)))
  integratedDF <- rbind(
    Beltran_NEvsAdeno_gExpr_sig_df_cisGenes,
    SU2CEastCoast_NEvsAdeno_gExpr_sig_df_cisGenes,
    SU2CWestCoast_NEvsAdeno_gExpr_sig_df_cisGenes
  )
  return(integratedDF)
}
mergeMatrices <- function(mat1, mat2){
  shared_rownames <- intersect(rownames(mat1), rownames(mat2))
  mat1 <- mat1[shared_rownames,]
  mat2 <- mat2[shared_rownames,]
  stopifnot(identical(rownames(mat1), rownames(mat2)))
  matM <- cbind(mat1, mat2)
  return(matM)
}
mean_of_3 <- function(x, y, z){
  mean(c(x, y, z))
}
integrateSigDFList <- function(SU2CEastCoast_NEvsAdeno_gExpr_sig_df_list){
  SU2CEastCoast_NEvsAdeno_gExpr_sig_df <- SU2CEastCoast_NEvsAdeno_gExpr_sig_df_list[[1]]
  SU2CEastCoast_NEvsAdeno_gExpr_sig_df$logFC <- NA
  SU2CEastCoast_NEvsAdeno_gExpr_sig_df$AveExpr <- NA
  SU2CEastCoast_NEvsAdeno_gExpr_sig_df$t <- NA
  SU2CEastCoast_NEvsAdeno_gExpr_sig_df$P.Value <- NA
  SU2CEastCoast_NEvsAdeno_gExpr_sig_df$adj.P.Val <- NA
  SU2CEastCoast_NEvsAdeno_gExpr_sig_df$B <- NA
  
  n_iters <- length(SU2CEastCoast_NEvsAdeno_gExpr_sig_df_list)
  
  my_colnames <- c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")
  for(i in 1:length(my_colnames)){
    col_i <- my_colnames[i]
    col_i_mat <- matrix(
      data = c(
        SU2CEastCoast_NEvsAdeno_gExpr_sig_df_list[[1]][[col_i]],
        SU2CEastCoast_NEvsAdeno_gExpr_sig_df_list[[2]][[col_i]],
        SU2CEastCoast_NEvsAdeno_gExpr_sig_df_list[[3]][[col_i]],
        SU2CEastCoast_NEvsAdeno_gExpr_sig_df_list[[4]][[col_i]],
        if(n_iters >= 5) {SU2CEastCoast_NEvsAdeno_gExpr_sig_df_list[[5]][[col_i]]} else {NULL},
        if(n_iters >= 6) {SU2CEastCoast_NEvsAdeno_gExpr_sig_df_list[[6]][[col_i]]} else {NULL},
        if(n_iters >= 7) {SU2CEastCoast_NEvsAdeno_gExpr_sig_df_list[[7]][[col_i]]} else {NULL},
        if(n_iters >= 8) {SU2CEastCoast_NEvsAdeno_gExpr_sig_df_list[[8]][[col_i]]} else {NULL},
        if(n_iters >= 9) {SU2CEastCoast_NEvsAdeno_gExpr_sig_df_list[[9]][[col_i]]} else {NULL},
        if(n_iters >= 10) {SU2CEastCoast_NEvsAdeno_gExpr_sig_df_list[[10]][[col_i]]} else {NULL}
      ),
      nrow = nrow(SU2CEastCoast_NEvsAdeno_gExpr_sig_df)
    )
    #if(col_i %in%  c("logFC", "AveExpr", "t", "B")){
    SU2CEastCoast_NEvsAdeno_gExpr_sig_df[[col_i]] <- rowMeans(col_i_mat)
    #} else {
    #for(j in 1:nrow(SU2CEastCoast_NEvsAdeno_gExpr_sig_df)){
    # SU2CEastCoast_NEvsAdeno_gExpr_sig_df[j, col_i] <- fishers_method(col_i_mat[j, ])
    #}
    #}
  }
  return(SU2CEastCoast_NEvsAdeno_gExpr_sig_df)
}

# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
# -----------------------------------------------------------------------------------
# -------------------------------- LOAD THE DATASETS --------------------------------
# -----------------------------------------------------------------------------------
# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&

Beltran_gExpr_tpm <- readRDS("data/Beltran-2016/beltran-sym-tpm-original.rds")
Beltran_metadata <- readRDS("experiments/oncomatch-analysis/processed_data/metadata/Beltran_metadata_forOMPlot.rds")
Beltran_metadata_NE_feats <- Beltran_metadata$Neuroendocrine_Features
names(Beltran_metadata_NE_feats) <- Beltran_metadata$Sample_Identifier
# Beltran_metadata_pathology <- Beltran_metadata$Pathology_Classification
# names(Beltran_metadata_pathology) <- Beltran_metadata$Sample_Identifier

SU2CEastCoast_gExpr_tpm <- readRDS("experiments/oncomatch-analysis/processed_data/gexpr/SU2CEastCoast_gExpr_polyA_tpm.rds")
SU2CEastCoast_metadata <- readRDS("experiments/oncomatch-analysis/processed_data/metadata/SU2C_EastCoast_metadata_forOMPlot.rds")
SU2CEastCoast_metadata <- SU2CEastCoast_metadata %>%
  filter(!is.na(Neuroendocrine_Features)) #%>%
  # filter(Pathology_Classification %in% c("Adenocarcinoma", "Pure small cell"))
SU2CEastCoast_gExpr_tpm <- SU2CEastCoast_gExpr_tpm[, SU2CEastCoast_metadata$Sample_Identifier]
SU2CEastCoast_metadata_NE_feats <- SU2CEastCoast_metadata$Neuroendocrine_Features
names(SU2CEastCoast_metadata_NE_feats) <- SU2CEastCoast_metadata$Sample_Identifier
# SU2CEastCoast_metadata_pathology <- SU2CEastCoast_metadata$Pathology_Classification
# names(SU2CEastCoast_metadata_pathology) <- SU2CEastCoast_metadata$Sample_Identifier

SU2CWestCoast_gExpr_tpm <- readRDS("data/SU2C_WestCoast/SU2CWestCoast_merged_tpm_geneNames.rds")
colnames(SU2CWestCoast_gExpr_tpm) = colnames(SU2CWestCoast_gExpr_tpm) %>%
  stringr::str_replace_all("_2018", "") %>%
  stringr::str_replace_all("\\.", "-")
SU2CWestCoast_gExpr_counts <- readRDS("data/SU2C_WestCoast/SU2CWestCoast_merged_counts.rds")
SU2CWestCoast_metadata <- readRDS("experiments/oncomatch-analysis/processed_data/metadata/SU2C_WestCoast_metadata_forOMPlot.rds")
SU2CWestCoast_metadata <- SU2CWestCoast_metadata %>%
  filter(!is.na(Neuroendocrine_Features))# %>%
  # filter(Pathology_Classification %in% c("Adenocarcinoma", "Pure small cell"))
SU2CWestCoast_gExpr_tpm <- SU2CWestCoast_gExpr_tpm[, SU2CWestCoast_metadata$Sample_Identifier]
SU2CWestCoast_metadata_NE_feats <- SU2CWestCoast_metadata$Neuroendocrine_Features
names(SU2CWestCoast_metadata_NE_feats) <- SU2CWestCoast_metadata$Sample_Identifier
# SU2CWestCoast_metadata_pathology <- SU2CWestCoast_metadata$Pathology_Classification
# names(SU2CWestCoast_metadata_pathology) <- SU2CWestCoast_metadata$Sample_Identifier

# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
# -----------------------------------------------------------------------------------
# ----------------------- COMPUTE THE DIFFERENTIAL EXPRESSION -----------------------
# -----------------------------------------------------------------------------------
# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
Beltran_NEvsAdeno_gExpr_sig_df <- runLimma(gExpr_tpm = Beltran_gExpr_tpm,
                                           NE_feats = Beltran_metadata_NE_feats)

SU2CEastCoast_NEvsAdeno_gExpr_sig_df_list <- list()

SU2CEastCoast_metadata_no_NE_feats_sampleNames <-
  names(SU2CEastCoast_metadata_NE_feats)[SU2CEastCoast_metadata_NE_feats=="No"]
SU2CEastCoast_metadata_no_NE_feats_sampleNames_list <- split(SU2CEastCoast_metadata_no_NE_feats_sampleNames, 1:10)
SU2CEastCoast_metadata_yes_NE_feats_sampleNames <-
  names(SU2CEastCoast_metadata_NE_feats)[SU2CEastCoast_metadata_NE_feats=="Yes"]
for(i in 1:10){
  # 210/21=10
  # SU2CEastCoast_NEvsAdeno_gExpr_sig_df <- runLimma(gExpr_tpm = SU2CEastCoast_gExpr_tpm,
  #                                                  NE_feats = SU2CEastCoast_metadata_NE_feats) %>%
  #   dplyr::rename("ID" = "rowname")
  this_iter_names <- c(SU2CEastCoast_metadata_yes_NE_feats_sampleNames,
                       SU2CEastCoast_metadata_no_NE_feats_sampleNames_list[[i]])
  SU2CEastCoast_NEvsAdeno_gExpr_sig_df_list[[i]] <- runLimma(gExpr_tpm = SU2CEastCoast_gExpr_tpm[, this_iter_names],
                                                             NE_feats = SU2CEastCoast_metadata_NE_feats[this_iter_names]) %>%
    dplyr::rename("ID" = "rowname")
  # Have all the tibbles have the same order as #1
  SU2CEastCoast_NEvsAdeno_gExpr_sig_df_list[[i]] <- SU2CEastCoast_NEvsAdeno_gExpr_sig_df_list[[i]] %>% arrange(match(SU2CEastCoast_NEvsAdeno_gExpr_sig_df_list[[i]]$ID, SU2CEastCoast_NEvsAdeno_gExpr_sig_df_list[[1]]$ID))
}

SU2CEastCoast_NEvsAdeno_gExpr_sig_df <- integrateSigDFList(SU2CEastCoast_NEvsAdeno_gExpr_sig_df_list)





SU2CWestCoast_NEvsAdeno_gExpr_sig_df_list <- list()
SU2CWestCoast_metadata_no_NE_feats_sampleNames <-
  names(SU2CWestCoast_metadata_NE_feats)[SU2CWestCoast_metadata_NE_feats=="No"]
SU2CWestCoast_metadata_no_NE_feats_sampleNames_list <- split(SU2CWestCoast_metadata_no_NE_feats_sampleNames, 1:4)
SU2CWestCoast_metadata_yes_NE_feats_sampleNames <-
  names(SU2CWestCoast_metadata_NE_feats)[SU2CWestCoast_metadata_NE_feats=="Yes"]

for(i in 1:4){
  # 68/17=4
  # SU2CWestCoast_NEvsAdeno_gExpr_sig_df <- runLimma(gExpr_tpm = SU2CWestCoast_gExpr_tpm,
  #                                                  NE_feats = SU2CWestCoast_metadata_NE_feats) %>%
  #   dplyr::rename("ID" = "rowname")
  this_iter_names <- c(SU2CWestCoast_metadata_yes_NE_feats_sampleNames,
                       SU2CWestCoast_metadata_no_NE_feats_sampleNames_list[[i]])
  SU2CWestCoast_NEvsAdeno_gExpr_sig_df_list[[i]] <- runLimma(gExpr_tpm = SU2CWestCoast_gExpr_tpm[, this_iter_names],
                                                             NE_feats = SU2CWestCoast_metadata_NE_feats[this_iter_names]) %>%
    dplyr::rename("ID" = "rowname")
  SU2CWestCoast_NEvsAdeno_gExpr_sig_df_list[[i]] <- SU2CWestCoast_NEvsAdeno_gExpr_sig_df_list[[i]] %>% arrange(match(SU2CWestCoast_NEvsAdeno_gExpr_sig_df_list[[i]]$ID, SU2CWestCoast_NEvsAdeno_gExpr_sig_df_list[[1]]$ID))
}
SU2CWestCoast_NEvsAdeno_gExpr_sig_df <- integrateSigDFList(SU2CWestCoast_NEvsAdeno_gExpr_sig_df_list)

cis_genes_top <- getTopCISGenes()
cis_genes_all <- getAllCISGenes()

Beltran_NEvsAdeno_gExpr_sig_df_cisGenes <- Beltran_NEvsAdeno_gExpr_sig_df %>%
  filter(ID %in% cis_genes_all)
SU2CEastCoast_NEvsAdeno_gExpr_sig_df_cisGenes <- SU2CEastCoast_NEvsAdeno_gExpr_sig_df %>%
  filter(ID %in% cis_genes_all)
SU2CWestCoast_NEvsAdeno_gExpr_sig_df_cisGenes <- SU2CWestCoast_NEvsAdeno_gExpr_sig_df %>%
  filter(ID %in% cis_genes_all)






combined_gExpr_sig_df_cisGenes_All <- 
  dplyr::left_join(
    rename(select(Beltran_NEvsAdeno_gExpr_sig_df_cisGenes, "ID", "adj.P.Val"), "adj.P.Val_Beltran" = "adj.P.Val"),
    rename(select(SU2CEastCoast_NEvsAdeno_gExpr_sig_df_cisGenes, "ID", "adj.P.Val"), "adj.P.Val_SU2C_EC" = "adj.P.Val"),
    by = c("ID" = "ID")
  ) %>%
  dplyr::left_join(
    rename(select(SU2CWestCoast_NEvsAdeno_gExpr_sig_df_cisGenes, "ID", "adj.P.Val"), "adj.P.Val_SU2C_WC" = "adj.P.Val"),
    by = c("ID" = "ID")
  )
for(i in 1:nrow(combined_gExpr_sig_df_cisGenes_All)){
  combined_gExpr_sig_df_cisGenes_All$fischer_integrated_adj.P.Val[i] <-
    fishers_method_for_mutate(
      p1 = combined_gExpr_sig_df_cisGenes_All$adj.P.Val_Beltran[i],
      p2 = combined_gExpr_sig_df_cisGenes_All$adj.P.Val_SU2C_EC[i],
      p3 = combined_gExpr_sig_df_cisGenes_All$adj.P.Val_SU2C_WC[i])
}

combined_gExpr_sig_df_cisGenes_topHits <- combined_gExpr_sig_df_cisGenes_All %>%
  filter(ID %in% cis_genes_top)

integratedDF_cisGenes_All <- getIntegratedGExprSigDFCISGenes(
  Beltran_NEvsAdeno_gExpr_sig_df_cisGenes,
  SU2CEastCoast_NEvsAdeno_gExpr_sig_df_cisGenes,
  SU2CWestCoast_NEvsAdeno_gExpr_sig_df_cisGenes,
  combined_gExpr_sig_df_cisGenes_All
)

integratedDF_cisGenes_Top <- integratedDF_cisGenes_All %>%
  filter(ID %in% cis_genes_top)

Beltran_NEvsAdeno_gExpr_sig <- matrix(Beltran_NEvsAdeno_gExpr_sig_df$AveExpr, ncol = 1)
rownames(Beltran_NEvsAdeno_gExpr_sig) <- Beltran_NEvsAdeno_gExpr_sig_df$ID
SU2CEastCoast_NEvsAdeno_gExpr_sig <- matrix(SU2CEastCoast_NEvsAdeno_gExpr_sig_df$AveExpr, ncol = 1)
rownames(SU2CEastCoast_NEvsAdeno_gExpr_sig) <- SU2CEastCoast_NEvsAdeno_gExpr_sig_df$ID
SU2CWestCoast_NEvsAdeno_gExpr_sig <- matrix(SU2CWestCoast_NEvsAdeno_gExpr_sig_df$AveExpr, ncol = 1)
rownames(SU2CWestCoast_NEvsAdeno_gExpr_sig) <- SU2CWestCoast_NEvsAdeno_gExpr_sig_df$ID

su2cEC_regulon <- readRDS("data/SU2C_EastCoast/SU2CEastCoast_net_pruned.rds")
su2cWC_regulon <- readRDS("data/SU2C_WestCoast/SU2CWestCoast_net_pruned.rds")

Beltran_NEvsAdeno_pAct_sig <- viper::viper(
  eset = Beltran_NEvsAdeno_gExpr_sig,
  regulon = list(su2cEC_regulon, su2cWC_regulon),
  method = 'none',
  eset.filter = FALSE,
  mvws = 10
)
SU2CEastCoast_NEvsAdeno_pAct_sig <- viper::viper(
  eset = SU2CEastCoast_NEvsAdeno_gExpr_sig,
  regulon = list(su2cEC_regulon, su2cWC_regulon),
  method = 'none',
  eset.filter = FALSE,
  mvws = 10
)
SU2CWestCoast_NEvsAdeno_pAct_sig <- viper::viper(
  eset = SU2CWestCoast_NEvsAdeno_gExpr_sig,
  regulon = list(su2cEC_regulon, su2cWC_regulon),
  method = 'none',
  eset.filter = FALSE,
  mvws = 10
)


Combined_NEvsAdeno_pAct_sig <- data.frame(
  "Beltran_pAct" = as.numeric(Beltran_NEvsAdeno_pAct_sig),
  "SU2CEastCoast_pAct" = SU2CEastCoast_NEvsAdeno_pAct_sig[rownames(Beltran_NEvsAdeno_pAct_sig), ],
  "SU2CWestCoast_pAct" = SU2CWestCoast_NEvsAdeno_pAct_sig[rownames(Beltran_NEvsAdeno_pAct_sig), ]
) %>%
  rownames_to_column("ID") %>%
  as_tibble()
for(i in 1:nrow(Combined_NEvsAdeno_pAct_sig)){
  Combined_NEvsAdeno_pAct_sig$mean_pAct[i] <-
    mean_of_3(
      Combined_NEvsAdeno_pAct_sig$Beltran_pAct[i],
      Combined_NEvsAdeno_pAct_sig$SU2CEastCoast_pAct[i],
      Combined_NEvsAdeno_pAct_sig$SU2CWestCoast_pAct[i])
}

Combined_NEvsAdeno_pAct_sig_cisGenes_All <- Combined_NEvsAdeno_pAct_sig %>%
  filter(ID %in% cis_genes_all) %>%
  arrange(desc(Beltran_pAct))
Combined_NEvsAdeno_pAct_sig_cisGenes_Top <- Combined_NEvsAdeno_pAct_sig %>%
  filter(ID %in% cis_genes_top) %>%
  arrange(desc(Beltran_pAct))


# View(combined_gExpr_sig_df_cisGenes_All)
# View(combined_gExpr_sig_df_cisGenes_topHits)
# View(Combined_NEvsAdeno_pAct_sig_cisGenes_All)
# View(Combined_NEvsAdeno_pAct_sig_cisGenes_Top)
# View(integratedDF_cisGenes_All)
# View(integratedDF_cisGenes_Top)

write.csv(arrange(combined_gExpr_sig_df_cisGenes_All, fischer_integrated_adj.P.Val), "experiments/NEPC-sig-analysis/reports/combined_gExpr_sig_df_cisGenes_All.csv")
write.csv(arrange(combined_gExpr_sig_df_cisGenes_topHits, fischer_integrated_adj.P.Val), "experiments/NEPC-sig-analysis/reports/combined_gExpr_sig_df_cisGenes_topHits.csv")

write.csv(arrange(Combined_NEvsAdeno_pAct_sig_cisGenes_All, mean_pAct), "experiments/NEPC-sig-analysis/reports/Combined_NEvsAdeno_pAct_sig_cisGenes_All.csv")
write.csv(arrange(Combined_NEvsAdeno_pAct_sig_cisGenes_Top, mean_pAct), "experiments/NEPC-sig-analysis/reports/Combined_NEvsAdeno_pAct_sig_cisGenes_Top.csv")

write.csv(arrange(integratedDF_cisGenes_Top, fischer_integrated_adj.P.Val), "experiments/NEPC-sig-analysis/reports/integratedDF_cisGenes_Top.csv")
write.csv(arrange(integratedDF_cisGenes_All, fischer_integrated_adj.P.Val), "experiments/NEPC-sig-analysis/reports/integratedDF_cisGenes_All.csv")

