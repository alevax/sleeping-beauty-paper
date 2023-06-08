setwd("~/Desktop/Califano_Lab_Fall_2021/sleeping-beauty-paper/")

# Libraries
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(magick))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(viridisLite))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(viper))

# Functions
source("libs/tools/oncomatch_funcs.R")
# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
# -----------------------------------------------------------------------------------
# ------------------------------- MERGE MATRICES FUNC -------------------------------
# -----------------------------------------------------------------------------------
# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
mergeMatrices <- function(mat1, mat2){
  shared_rownames <- intersect(rownames(mat1), rownames(mat2))
  mat1 <- mat1[shared_rownames,]
  mat2 <- mat2[shared_rownames,]
  matM <- cbind(mat1, mat2)
  print(table(colnames(matM) == c(colnames(mat1), colnames(mat2))))
  colnames(matM) <- c(colnames(mat1), colnames(mat2))
  return(matM)
}
range_standardize <- function(x, newMin = 0, newMax = 1){
  x <- log(x + 1, base = 2)
  (x - min(x))/(max(x)-min(x)) * (newMax - newMin) + newMin 
}
# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
# -----------------------------------------------------------------------------------
# ------------------------------ ADD TO METADATA FUNCS ------------------------------
# -----------------------------------------------------------------------------------
# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
addGeneToMetadataDF <- function(metadata_df, gene_name, gExpr){
  if(gene_name %in% rownames(gExpr)){
    metadata_df[[paste0(gene_name, "_gExpr")]] <-
      as.numeric(gExpr[gene_name, metadata_df$Sample_Identifier])
  }
  return(metadata_df)
}
addProteinToMetadataDF <- function(metadata_df, protein_name, pAct){
  if(protein_name %in% rownames(pAct)){
    metadata_df[[paste0(protein_name, "_pAct")]] <-
      as.numeric(pAct[protein_name, metadata_df$Sample_Identifier])
  }
  return(metadata_df)
}
# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
# -----------------------------------------------------------------------------------
# ---------------- HIEARCHICAL CLUSTERING OF OM MATRIX ROWS & COLUMNS ---------------
# -----------------------------------------------------------------------------------
# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&

# --------------------------------- HELPER FUNCTION ---------------------------------
HClustPatientOMDataset <- function(Human_om, Human_merge_metadata, dataset, clust_method = "ward.D2"){
  Human_om_dataset <- Human_om[,Human_merge_metadata$Dataset == dataset]
  set.seed(0)
  hclust_order <- hclust(d = as.dist(1-cor(Human_om_dataset)),method = clust_method)$order
  return(colnames(Human_om_dataset)[hclust_order])
}

# ---------------------------------- MAIN FUNCTIONS ---------------------------------
# Sort the OM rows by using hierarchical clustering
clusterOMRowsByHClustering <- function(Human_om, gemm_metadata, clust_method = "ward.D2"){
  # Ensure the order is the same
  gemm_metadata <- gemm_metadata %>%
    arrange(match(Sample_Identifier, rownames(Human_om))) 
  # HClust
  set.seed(0)
  hclust_order <- hclust(d = as.dist(1-cor(t(Human_om))), method = clust_method)$order
  samplesOrdered <- rownames(Human_om)[hclust_order]
  return(Human_om[samplesOrdered,])
}
# Sort the OM columns by using hierarchical clustering
clusterOMColsByHClustering <- function(Human_om, Human_merge_metadata, clust_method = "ward.D2"){
  # Ensure the order is the same
  Human_merge_metadata <- Human_merge_metadata %>%
    dplyr::arrange(match(Sample_Identifier, colnames(Human_om))) 
  # HClust each group
  EC_hclust_samplesOrdered <- HClustPatientOMDataset(Human_om,
                                                     Human_merge_metadata,
                                                     dataset = "SU2C_East_Coast",
                                                     clust_method = clust_method)
  WC_hclust_samplesOrdered <- HClustPatientOMDataset(Human_om,
                                                     Human_merge_metadata,
                                                     dataset = "SU2C_West_Coast",
                                                     clust_method = clust_method)
  Beltran_hclust_samplesOrdered <- HClustPatientOMDataset(Human_om,
                                                          Human_merge_metadata,
                                                          dataset = "Beltran",
                                                          clust_method = clust_method)
  samplesOrdered <- c(
    EC_hclust_samplesOrdered,
    WC_hclust_samplesOrdered,
    Beltran_hclust_samplesOrdered
  )
  return(Human_om[,samplesOrdered])
}

# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
# -----------------------------------------------------------------------------------
# ------------------------------ GET HEATMAP LIST FUNCS -----------------------------
# -----------------------------------------------------------------------------------
# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&

# --------------------------------- HELPER FUNCTION ---------------------------------
getVega20 <- function(){
  vega_10_dark = c('#1f77b4',
                   '#ff7f0e',
                   '#2ca02c',
                   '#d62728',
                   '#9467bd',
                   '#8c564b',
                   '#e377c2',
                   '#7f7f7f',
                   '#bcbd22',
                   '#17becf')
  vega_10_light = c('#aec7e8',
                    '#ffbb78',
                    '#98df8a',
                    '#ff9896',
                    '#c5b0d5',
                    '#c49c94',
                    '#f7b6d2',
                    '#c7c7c7',
                    '#dbdb8d',
                    '#9edae5')
  vega_20_complete = c(vega_10_dark, vega_10_light)
  vega_20_complete
}
getNPats2ClonesDFAnnotRowsRight <- function(om, om_min_threshold){
  # Create the n_pats2clones annotation
  n_pats2clones <- rowSums(om > om_min_threshold)
  df_annot_rows_right <- rowAnnotation( n_pats2clones = anno_barplot(
    n_pats2clones,
    bar_width = 1,
    gp = gpar(col = "white", fill = "#FFE200"),
    border = FALSE,
    width = unit(2, "cm") )
  )
  return(df_annot_rows_right)
}
getNClones2PatDFAnnotColsTop <- function(om, om_min_threshold){
  # Create the n_clones2pat annotation
  n_clones2pat <- colSums(om > om_min_threshold)
  df_annot_cols_top <-  HeatmapAnnotation( n_clones2pat = anno_barplot(
    n_clones2pat ,
    bar_width = 1,
    gp = gpar(col = "white", fill = "#FFE200"),
    border = FALSE,
    width = unit(2, "cm") )
  )
  return(df_annot_cols_top)
}
getDfAnnotRowsLeft <- function(df_annot_rows){
  phenotype_color.vector <- rev(getVega20()[1:2])
  names(phenotype_color.vector) <- names(table(df_annot_rows$Phenotype))
  
  sb_status_color.vector <- c("red", "white")
  names(sb_status_color.vector) <- names(table(df_annot_rows$SB_Status))
  
  df_annot_rows_left <- 
    rowAnnotation(
      df = df_annot_rows, #CANNOT HAVE THIS AS A TIBBLE
      col = list(
        Phenotype = phenotype_color.vector,
        SB_Status = sb_status_color.vector
      ),
      annotation_legend_param = list(
        Phenotype = list(
          title_gp = gpar(fontsize = 20, fontface = "bold"),
          labels_gp = gpar(fontsize = 20)
        ),
        SB_Status = list(
          title_gp = gpar(fontsize = 20, fontface = "bold"),
          labels_gp = gpar(fontsize = 20)
        )
      ),
      annotation_name_gp= gpar(fontsize = 20),
      simple_anno_size = unit(x = 10, units = "mm")
    )
  return(df_annot_rows_left)
}
getDfAnnotCols <- function(df_annot_filtered, extra_annots = TRUE, pathology = FALSE, annot_fontsize = 20){
  require(circlize)
  require(viridis)
  my_col_list <- list(
    Dataset = c(
      "SU2C_East_Coast" = "#003884",
      "SU2C_West_Coast" = "#B71234",
      "Beltran" = "#E09F3E"# "goldenrod3" #"#FFD700"
    ),
    Tissue_Site = c(
      "Adrenal" = "red",
      "Bone" = "black" ,
      "LN" = "orange" ,
      "Other Soft tissue" = "white",
      "Liver" = "violet",
      "Bladder" = "green",
      "Lung" = "blue",
      "Prostate" = "yellowgreen",
      "Brain" = "pink",
      "Skull base" = "darkgray",
      "Pelvic mass" = "coral"
    ),
    Pathology_Classification = c(
      "Large cell NE carcinoma" = "deeppink1",
      "Pure small cell" = "red",
      "Mixed small cell" = "darkgoldenrod1",
      "Not small cell" =  '#63666a',
      "Unclassified" = "lightgray",
      "Adenocarcinoma with NE features" = "darkorange4",
      "Adenocarcinoma" =  '#63666a'
    ),
    Neuroendocrine_Features = c("No" = "white",
                                "Yes" = '#1f77b4'),#'red'),
    NEPC_Score = colorRamp2(seq(min(df_annot_filtered$NEPC_Score, na.rm = TRUE),
                                max(df_annot_filtered$NEPC_Score, na.rm = TRUE),
                                length.out = 10) , inferno(10)),
    KLK3_gExpr = colorRamp2(seq(
      0,
      summary(df_annot_filtered$KLK3_gExpr)["3rd Qu."],
      length.out = 10),
      viridis(10)),
    SYP_gExpr = colorRamp2(seq(0,
                               max(df_annot_filtered$SYP_gExpr),
                               length.out = 10),
                           viridis(10)),
    CHGA_gExpr = colorRamp2(seq(0,
                                max(df_annot_filtered$CHGA_gExpr),
                                length.out = 10),
                            viridis(10)),
    AR_pAct = colorRamp2(c(-10, 0, 5), c("blue", "white", "red")),
    NSD2_pAct = colorRamp2(c(-10, 0, 5), c("blue", "white", "red"))
  )
  my_annotation_legend_params = list(
    Dataset = list(
      title_gp = gpar(fontsize = annot_fontsize, fontface = "bold"),
      labels_gp = gpar(fontsize = annot_fontsize)
    ),
    Tissue_Site = list(
      title_gp = gpar(fontsize = annot_fontsize, fontface = "bold"),
      labels_gp = gpar(fontsize = annot_fontsize)
    ),
    Neuroendocrine_Features = list(
      title_gp = gpar(fontsize = annot_fontsize, fontface = "bold"),
      labels_gp = gpar(fontsize = annot_fontsize)
    ),
    NEPC_Score = list(
      title_gp = gpar(fontsize = annot_fontsize, fontface = "bold"),
      labels_gp = gpar(fontsize = annot_fontsize)
    )
  )
  df_annot_cols = HeatmapAnnotation(
    df = data.frame(dplyr::select(df_annot_filtered,
                                  "Dataset",
                                  "Tissue_Site",
                                  "Neuroendocrine_Features",
                                  "NEPC_Score",
                                  if(pathology) {"Pathology_Classification"} else {NULL},
                                  if(extra_annots) {"KLK3_gExpr"} else {NULL},
                                  if(extra_annots) {"SYP_gExpr"} else {NULL},
                                  if(extra_annots) {"CHGA_gExpr"} else {NULL},
                                  if(extra_annots) {"AR_pAct"} else {NULL},
    )),
    col = my_col_list,
    annotation_name_gp= gpar(fontsize = annot_fontsize),
    simple_anno_size = unit(x = 10, units = "mm"),
    annotation_legend_param = my_annotation_legend_params,
    show_legend = c(TRUE,
                    TRUE,
                    TRUE,
                    TRUE,
                    if(pathology) {FALSE} else {NULL},
                    if(extra_annots) {FALSE} else {NULL},
                    if(extra_annots) {FALSE} else {NULL},
                    if(extra_annots) {FALSE} else {NULL},
                    if(extra_annots) {FALSE} else {NULL}
    )
  )
  return(df_annot_cols)
}
getMainHeatmap <- function(om,
                           df_annot_rows_left,
                           df_annot_rows_right,
                           Human_merge_col_split,
                           df_annot_rows){
  om_binned <- OncoLoopScoresBinning(om, cutoff_max = 35)
  main_hm_binned <- Heatmap(om_binned$mat,
                            col = om_binned$colors,
                            cluster_rows =  FALSE,
                            use_raster = TRUE ,
                            name = "OncoMatch",
                            rect_gp = gpar(col = "white", lwd = 2),
                            left_annotation = df_annot_rows_left,
                            right_annotation = df_annot_rows_right ,
                            show_row_names = TRUE,
                            show_column_names = FALSE ,
                            row_names_gp = gpar(fontsize = 5),
                            cluster_row_slices = FALSE,
                            row_split = NULL,
                            row_gap = unit(0.15,"in"),
                            column_split = Human_merge_col_split,
                            column_gap = unit(0.9,"in"),
                            column_title_gp = gpar(fontsize = 25,
                                                   fontface = "bold",
                                                   col = c("#B71234", "#003884", "#E09F3E")),
                            border_gp = gpar(col = "lightgray", lty = 1, lwd = 3, alpha = 0.6),
                            heatmap_legend_param = list(
                              title_gp = gpar(fontsize = 20, fontface = "bold"),
                              labels_gp = gpar(fontsize = 20)
                            )
  )
  return(main_hm_binned)
}

# ---------------------------------- MAIN FUNCTION ----------------------------------
getHeatmapList <- function(Human_merge_om, gemm_metadata, Human_merge_metadata, om_min_threshold = 10){
  # Ensure the sample order in the PlateSeq Metadata match the row order of the OncoMatch Matrices
  gemm_df_annot_rows_Human_merge <- gemm_metadata %>%
    arrange(match(Sample_Identifier, rownames(Human_merge_om))) %>%
    dplyr::select(-c("Sample_Identifier")) %>%
    dplyr::select(
      "Phenotype",
      "SB_Status"
      ) %>%
    data.frame()
  
  # Use the ordered PlateSeq Metadata to compute the row annotations (on left of main heatmap)
  gemm_df_annot_rows_left_Human_merge <- getDfAnnotRowsLeft(df_annot_rows = gemm_df_annot_rows_Human_merge)
  
  # Ensure the sample order in the SU2C Metadata match the column order of the OncoMatch Matrices
  Human_merge_metadata <- Human_merge_metadata %>%
    arrange(match(Sample_Identifier, colnames(Human_merge_om)))
  
  # Use the SU2C Metadata to compute the column annotations
  Human_merge_df_annot_cols <- getDfAnnotCols(df_annot_filtered = Human_merge_metadata)
  Human_merge_df_annot_cols_simpleAnnot <- getDfAnnotCols(df_annot_filtered = Human_merge_metadata, extra_annots = FALSE)
  
  # Use the OncoMatch Matrices to compute the # of patients & clones matching annotations
  gemm_df_annot_rows_Human_merge_right <-
    getNPats2ClonesDFAnnotRowsRight(Human_merge_om, om_min_threshold)
  Human_merge_df_annot_cols_top <-
    getNClones2PatDFAnnotColsTop(Human_merge_om, om_min_threshold)
  
  # Compute the main heatmaps
  Human_merge_col_split <- factor(Human_merge_metadata$Dataset,
                                  levels = c("Beltran",
                                             "SU2C_East_Coast",
                                             "SU2C_West_Coast"),
                                  ordered = TRUE)
  set.seed(0)
  Human_merge_main_hm <- getMainHeatmap(Human_merge_om,
                                        gemm_df_annot_rows_left_Human_merge,
                                        gemm_df_annot_rows_Human_merge_right,
                                        Human_merge_col_split,
                                        df_annot_rows)
  ht_list = Human_merge_df_annot_cols %v% Human_merge_df_annot_cols_top %v% Human_merge_main_hm
  ht_list_simpleAnnot = Human_merge_df_annot_cols_simpleAnnot %v% Human_merge_df_annot_cols_top %v% Human_merge_main_hm
  res <- list(ht_list, ht_list_simpleAnnot)
  names(res) <- c("ht_list", "ht_list_simpleAnnot")
  return(res)
}
drawHeatmapsLists <- function(res, out_dir){
  ht_opt(HEATMAP_LEGEND_PADDING = unit(1, "in"))
  p <- draw(res[["ht_list"]],
            column_km = 1,
            merge_legend = TRUE,
            heatmap_legend_side = "bottom", 
            annotation_legend_side = "bottom")#,
  # png(paste0(out_dir, "SleepingBeautyOncoMatch_extraAnnotations.png"), width = 32, height = 14 , unit="in", res=300)
  pdf(paste0(out_dir, "SleepingBeautyOncoMatch_extraAnnotations.pdf"), width = 70, height = 14)
  print(p)
  dev.off()
  p <- draw(res[["ht_list_simpleAnnot"]],
            column_km = 1,
            merge_legend = TRUE,
            heatmap_legend_side = "bottom", 
            annotation_legend_side = "bottom")#,
  # png(paste0(out_dir, "SleepingBeautyOncoMatch__simpleAnnot.png"), width = 32, height = 12 , unit="in", res=300)
  pdf(paste0(out_dir, "SleepingBeautyOncoMatch_simpleAnnot.pdf"), width = 70, height = 12)
  print(p)
  dev.off()
}

# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
# -----------------------------------------------------------------------------------
# ---------------------------- MAIN SAVE OM HEATMAP FUNC ----------------------------
# -----------------------------------------------------------------------------------
# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&

saveOMResultsWithThreshold <- function(Human_merge_om,
                                       gemm_metadata,
                                       Human_merge_metadata,
                                       out_dir,
                                       om_min_threshold = 0){
  Human_merge_om[Human_merge_om < om_min_threshold] <- 0
  
  # Add small almost-0 data for SU2C Samples with 0 values
  set.seed(0)
  index <- which(apply(Human_merge_om, 2, function(x)sum(x)==0))
  Human_merge_om[,index] <- runif(nrow(Human_merge_om),0,0.001)   
  
  set.seed(0)
  res_byHClust <- Human_merge_om %>%
    clusterOMRowsByHClustering(gemm_metadata, clust_method = "single") %>%
    clusterOMColsByHClustering(Human_merge_metadata, clust_method = "single") %>%
    getHeatmapList(gemm_metadata, Human_merge_metadata, om_min_threshold)
  drawHeatmapsLists(res_byHClust, out_dir = out_dir)
}

#####################################################################################
#####################################################################################
#################################                   #################################
#################################                   #################################
#################################                   ################################# 
################################# START OF ANALYSIS #################################
#################################                   #################################
#################################                   #################################
#################################                   #################################
#####################################################################################
#####################################################################################

# Load gExpr data
Su2c_EC_gExpr <- readRDS("experiments/oncomatch-analysis/processed_data/gexpr/SU2CEastCoast_gExpr_polyA_tpm.rds")
Su2c_WC_gExpr <- readRDS("data/SU2C_WestCoast/SU2CWestCoast_merged_tpm_geneNames.rds")
Beltran_gExpr <- readRDS("data/Beltran-2016/beltran-sym-tpm-original.rds")

# Load VIPER data
gemm_pAct <- readRDS("experiments/oncomatch-analysis/processed_data/viper/SleepingBeauty_viper_zscore_SU2CNets_SBNet.rds")
Su2c_EC_pAct <- readRDS("experiments/oncomatch-analysis/processed_data/viper/SU2CEastCoast_viper_zscore_SU2CNets.rds")
Su2c_WC_pAct <- readRDS("experiments/oncomatch-analysis/processed_data/viper/SU2CWestCoast_viper_zscore_SU2CNets.rds")
Beltran_pAct <- readRDS("experiments/oncomatch-analysis/processed_data/viper/Beltran_viper_zscore_SU2CNets.rds")

# Load OnCoMatch data
Su2c_EC_om <- readRDS("experiments/oncomatch-analysis/processed_data/oncomatch/SleepingBeauty_SU2CEastCoast_OncoMatch.rds")
Su2c_WC_om <- readRDS("experiments/oncomatch-analysis/processed_data/oncomatch/SleepingBeauty_SU2CWestCoast_OncoMatch.rds")
Beltran_om <- readRDS("experiments/oncomatch-analysis/processed_data/oncomatch/SleepingBeauty_Beltran_OncoMatch.rds")

# Load MetaData
Su2c_EC_metadata <- readRDS("experiments/oncomatch-analysis/processed_data/metadata/SU2C_EastCoast_metadata_forOMPlot.rds")
Su2c_WC_metadata <- readRDS("experiments/oncomatch-analysis/processed_data/metadata/SU2C_WestCoast_metadata_forOMPlot.rds")
Beltran_metadata <- readRDS("experiments/oncomatch-analysis/processed_data/metadata/Beltran_metadata_forOMPlot.rds")
gemm_metadata <- readRDS("experiments/oncomatch-analysis/processed_data/metadata/SleepingBeauty_metadata_forOMPlot.rds")


# Clean Up WC Sample Names For Plotting
colnames(Su2c_WC_gExpr) = colnames(Su2c_WC_gExpr) %>%
  stringr::str_replace_all("_2018", "") %>%
  stringr::str_replace_all("\\.", "-")
colnames(Su2c_WC_pAct) = colnames(Su2c_WC_pAct) %>%
  stringr::str_replace_all("_2018", "") %>%
  stringr::str_replace_all("\\.", "-")
colnames(Su2c_WC_om) = colnames(Su2c_WC_om) %>%
  stringr::str_replace_all("_2018", "") %>%
  stringr::str_replace_all("\\.", "-")


# Range standardize the gExpr
Su2c_EC_gExpr["KLK3",] <- range_standardize(Su2c_EC_gExpr["KLK3",])
Su2c_EC_gExpr["SYP",] <- range_standardize(Su2c_EC_gExpr["SYP",])
Su2c_EC_gExpr["CHGA",] <- range_standardize(Su2c_EC_gExpr["CHGA",])

Su2c_WC_gExpr["KLK3",] <- range_standardize(Su2c_WC_gExpr["KLK3",])
Su2c_WC_gExpr["SYP",] <- range_standardize(Su2c_WC_gExpr["SYP",])
Su2c_WC_gExpr["CHGA",] <- range_standardize(Su2c_WC_gExpr["CHGA",])

Beltran_gExpr["KLK3",] <- range_standardize(Beltran_gExpr["KLK3",])
Beltran_gExpr["SYP",] <- range_standardize(Beltran_gExpr["SYP",])
Beltran_gExpr["CHGA",] <- range_standardize(Beltran_gExpr["CHGA",])


# Merge data together
Human_merge_gExpr <- Su2c_EC_gExpr %>%
  mergeMatrices(Su2c_WC_gExpr) %>%
  mergeMatrices(Beltran_gExpr)
Human_merge_pAct <- Su2c_EC_pAct %>%
  mergeMatrices(Su2c_WC_pAct) %>%
  mergeMatrices(Beltran_pAct)
Human_merge_om <- Su2c_EC_om %>%
  mergeMatrices(Su2c_WC_om) %>%
  mergeMatrices(Beltran_om)
Human_merge_metadata <- Su2c_EC_metadata %>%
  dplyr::bind_rows(Su2c_WC_metadata) %>%
  dplyr::bind_rows(Beltran_metadata)
Human_merge_metadata$Dataset <- c(rep("SU2C_East_Coast", nrow(Su2c_EC_metadata)),
                                 rep("SU2C_West_Coast", nrow(Su2c_WC_metadata)),
                                 rep("Beltran", nrow(Beltran_metadata)))


# Add GExpr and PAct to metadata
genes_to_plot <- c("KLK3", "SYP", "CHGA")
proteins_to_plot <- c("AR", "NSD2")

for(i in 1:length(genes_to_plot)){
  Human_merge_metadata <- addGeneToMetadataDF(Human_merge_metadata,
                                              genes_to_plot[i],
                                              Human_merge_gExpr)
}
for(i in 1:length(proteins_to_plot)){
  Human_merge_metadata <- addProteinToMetadataDF(Human_merge_metadata,
                                                 proteins_to_plot[i],
                                                 Human_merge_pAct)
}


# Plot the Oncomatch Results
saveOMResultsWithThreshold(Human_merge_om,
                           gemm_metadata,
                           Human_merge_metadata,
                           out_dir = "experiments/oncomatch-analysis/reports/",
                           om_min_threshold = 5)


getTopCISGenes <- function(){
  top_hits_tbl <- readxl::read_xlsx("experiments/integrative-analysis/processed_data/top-hits-for-validation-final-fisher-integration-table.xlsx")
  cis_genes_top <- top_hits_tbl$gene_human
  return(cis_genes_top)
}
ZScoreTransform <- function(dat.mat, na.rm = TRUE){
  # generate GES
  ges.mat <- t(apply(dat.mat, 1, function(x) {
    (x - mean(x, na.rm = na.rm)) / sd(x, na.rm = na.rm)
  }))
  ges.mat
}
# source("~/Desktop/Califano_Lab_Fall_2021/sleeping-beauty-paper/libs/network_analysis.R")
getMyTable <- function(filename){
  my_table <- readRDS(filename)
  table(my_table$regulator_type)
  names(my_table)[names(my_table) == "gene_human"] <- "human_gene_name"
  return(my_table)
}
getTFsWithHighestNEPCScoresDF <- function(my_table, n_top = 50){
  # x = the TFs with the top 50 NEPC VIPER Scores
  top50_TF_df <- my_table %>% 
    filter( regulator_type == "TF" ) %>%
    arrange( desc(NEPC_viper_score) ) %>%
    slice_head(n = n_top)
  return(top50_TF_df)
}
getTFsWithLowestNEPCScoresDF <- function(my_table, n_lowest = 50){
  # y = the TFs with the lowest 50 NEPC VIPER Scores
  lowest50_TF_df <- my_table %>% 
    filter( regulator_type == "TF" ) %>%
    arrange( NEPC_viper_score ) %>%
    slice_head(n = n_lowest)
  return(lowest50_TF_df)
}
saveVIPERTopCISGenesHeatmap <- function(Su2c_EC_pAct,
                                        Su2c_WC_pAct,
                                        Beltran_pAct,
                                        Human_merge_metadata,
                                        protein_set_name = "cis",
                                        gExpr = FALSE,
                                        out_file,
                                        Su2c_EC_gExpr = NULL,
                                        Su2c_WC_gExpr = NULL,
                                        Beltran_gExpr = NULL){
  if(protein_set_name == "cis"){
    protein_set <- getTopCISGenes() %>%
      intersect(rownames(Beltran_pAct))
    my_row_split <- NULL
  } else if(protein_set_name == "mrs"){
    my_table <-
      getMyTable("experiments/integrative-analysis/processed_data/tibble-of-sb-cis-nepc-viper-cindy-integrative-analysis.rds")
    top_candidate_mrs <- getTFsWithHighestNEPCScoresDF(my_table, n_top = 50)$human_gene_name %>%
      intersect(rownames(Beltran_pAct))
    bottom_candidate_mrs <- getTFsWithLowestNEPCScoresDF(my_table, n_lowest = 50)$human_gene_name %>%
      intersect(rownames(Beltran_pAct))
    protein_set <- c(top_candidate_mrs, bottom_candidate_mrs)
    my_row_split <- factor(
      x = c(rep("top_mrs", length(top_candidate_mrs)), rep("bottom_mrs", length(bottom_candidate_mrs))),
      levels = c("top_mrs", "bottom_mrs"),
      ordered = TRUE
    )
  } else {
    stop("wrong protein set name")
  }
  
  Su2c_EC_pAct_CIS <- Su2c_EC_pAct[protein_set, ]
  Su2c_WC_pAct_CIS <- Su2c_WC_pAct[protein_set, ]
  Beltran_pAct_CIS <- Beltran_pAct[protein_set, ]
  
  if(gExpr){
    Su2c_EC_pAct_CIS <- ZScoreTransform(Su2c_EC_pAct_CIS) #t(scale(t(Su2c_EC_pAct_CIS)))
    Su2c_WC_pAct_CIS <- ZScoreTransform(Su2c_WC_pAct_CIS) #t(scale(t(Su2c_WC_pAct_CIS)))
    Beltran_pAct_CIS <- ZScoreTransform(Beltran_pAct_CIS) #t(scale(t(Beltran_pAct_CIS)))
  }
  
  Human_merge_pAct <- Su2c_EC_pAct_CIS %>%
    mergeMatrices(Su2c_WC_pAct_CIS) %>%
    mergeMatrices(Beltran_pAct_CIS)
  Human_merge_pAct <- Human_merge_pAct[, Human_merge_metadata$Sample_Identifier]
  
  Human_merge_pAct_hclust <- Human_merge_pAct %>%
    clusterOMColsByHClustering(Human_merge_metadata, clust_method = "single")
  Human_merge_metadata <- Human_merge_metadata %>%
    arrange(match(Sample_Identifier, colnames(Human_merge_om)))
  
  # Compute the main heatmaps
  
  NE_text_add_on <- Human_merge_metadata$Neuroendocrine_Features %>%
    stringr::str_replace_all("Yes", "_NEPC") %>%
    stringr::str_replace_all("No", "_CRPC")
  NE_text_add_on[is.na(NE_text_add_on)] <- "_CRPC"
  Human_merge_metadata$Dataset_NEPC <- factor(
    paste0(Human_merge_metadata$Dataset, NE_text_add_on),
    levels = c(
      "Beltran_CRPC",
      "Beltran_NEPC",
      "SU2C_East_Coast_CRPC",
      "SU2C_East_Coast_NEPC",
      "SU2C_West_Coast_CRPC",
      "SU2C_West_Coast_NEPC"
    ),
    ordered = TRUE
  )
  Human_merge_df_annot_cols_simpleAnnot <- getDfAnnotCols(df_annot_filtered = Human_merge_metadata,
                                                          extra_annots = FALSE,
                                                          pathology = FALSE,
                                                          annot_fontsize = ifelse(protein_set_name == "cis", 20, 35))
  
  my_split <- factor(Human_merge_metadata$Dataset_NEPC,
                     levels = c(
                       "Beltran_CRPC",
                       "Beltran_NEPC",
                       "SU2C_East_Coast_CRPC",
                       "SU2C_East_Coast_NEPC",
                       "SU2C_West_Coast_CRPC",
                       "SU2C_West_Coast_NEPC"
                     ),
                     ordered = TRUE
  )
  
  if(gExpr == TRUE){
    # my_colors <-  colorRamp2(seq(-4, 4, length = 3), c("purple", "#EEEEEE", "orange"))
    # my_colors <-  colorRamp2(seq(-4, 4, length = 5), c("magenta4", "magenta", "#EEEEEE", "springgreen", "springgreen4"))
    # my_colors <-  colorRamp2(seq(-4, 4, length = 3), c("magenta4", "#EEEEEE", "springgreen4"))
    # wistia <- c(
    #   "#e4ff7a",
    #   "#f1f34a",
    #   "#fbeb24",
    #   "#ffdd13",
    #   "#ffcc07",
    #   "#ffbb00",
    #   "#ffb100",
    #   "#ffa501",
    #   "#fe9a02",
    #   "#fd8e00",
    #   "#fc7f00"
    # )
    # my_colors <-  colorRamp2(seq(-4, 4, length = 11), wistia)
    # my_colors <-  colorRamp2(seq(-4, 4, length = 11), wistia)
    
    # PuOr <- c("#7F3B08", "#B35806", "#E08214", "#FDB863", "#FEE0B6", "#F7F7F7", "#D8DAEB", "#B2ABD2", "#8073AC", "#542788", "#2D004B")
    # PuOr2 <- c("#7F3B08", "#B35806", "#E08214", "#FDB863", "#FEE0B6", "#F7F7F7", "orchid", "orchid1", "orchid2", "orchid3", "orchid4")
    # PuOr3 <- c("#7F3B08", "#B35806", "#E08214", "#FDB863", "#FEE0B6", "#F7F7F7", "magenta", "magenta1", "magenta2", "magenta3", "magenta4")
    # my_colors <-  colorRamp2(seq(-4, 4, length = 7), rev(PuOr[3:9]))
    # my_colors <-  colorRamp2(seq(-4, 4, length = 7), rev(PuOr2[3:9]))
    my_colors <-  colorRamp2(seq(-4, 4, length = 16), viridis(16))
    # PiYG <- c("#8E0152", "#C51B7D", "#DE77AE", "#F1B6DA", "#FDE0EF", "#F7F7F7", "#E6F5D0", "#B8E186", "#7FBC41", "#4D9221", "#276419")
    # my_colors <- colorRamp2(seq(-4, 4, length = 11), colors = rev(PiYG))
    # my_colors <- NULL
    
  } else {
    my_colors <- NULL
  }
  
  set.seed(0)
  main_hm_pAct <- Heatmap(Human_merge_pAct_hclust,
                          col = my_colors,
                          cluster_columns =  TRUE,
                          cluster_rows =  TRUE,
                          use_raster = TRUE,
                          name = "Top CIS Hits VIPER",
                          rect_gp = gpar(col = "white", lwd = 2),
                          show_row_names = TRUE,
                          show_column_names = FALSE ,
                          row_names_gp = gpar(fontsize = ifelse(protein_set_name == "cis", 25, 45)),
                          cluster_row_slices = FALSE,
                          cluster_column_slices = FALSE,
                          row_split = my_row_split,
                          row_gap = unit(0.9,"in"),
                          column_split = my_split, #Human_merge_metadata$Dataset_NEPC,
                          column_gap = unit(0.9,"in"),
                          column_title_gp = gpar(fontsize = ifelse(protein_set_name == "cis", 25, 60),
                                                 fontface = "bold",
                                                 col = c("#B71234", "#B71234", "#003884", "#003884", "#E09F3E", "#E09F3E")),
                          border_gp = gpar(col = "lightgray", lty = 1, lwd = 3, alpha = 0.6),
                          heatmap_legend_param = list(
                            title_gp = gpar(fontsize = ifelse(protein_set_name == "cis", 20, 35), fontface = "bold"),
                            labels_gp = gpar(fontsize = ifelse(protein_set_name == "cis", 20, 35))
                          )
  )
  
  if(!is.null(Su2c_EC_gExpr)){
    Su2c_EC_gExpr_CIS <- ZScoreTransform(Su2c_EC_gExpr[protein_set, ])
    Su2c_WC_gExpr_CIS <- ZScoreTransform(Su2c_WC_gExpr[protein_set, ])
    Beltran_gExpr_CIS <- ZScoreTransform(Beltran_gExpr[protein_set, ])
    Human_merge_gExpr_CIS <- Su2c_EC_gExpr_CIS %>%
      mergeMatrices(Su2c_WC_gExpr_CIS) %>%
      mergeMatrices(Beltran_gExpr_CIS)
    Human_merge_gExpr_CIS <- Human_merge_gExpr_CIS[rownames(Human_merge_pAct_hclust), colnames(Human_merge_pAct_hclust)]
    set.seed(0)
    main_hm_gExpr <- Heatmap(Human_merge_gExpr_CIS,
                            col = colorRamp2(seq(-4, 4, length = 16), viridis(16)),
                            cluster_columns =  TRUE,
                            cluster_rows =  TRUE,
                            use_raster = TRUE,
                            name = "Top CIS Hits VIPER",
                            rect_gp = gpar(col = "white", lwd = 2),
                            show_row_names = TRUE,
                            show_column_names = FALSE ,
                            row_names_gp = gpar(fontsize = ifelse(protein_set_name == "cis", 25, 45)),
                            cluster_row_slices = FALSE,
                            cluster_column_slices = FALSE,
                            row_split = my_row_split,
                            row_gap = unit(0.9,"in"),
                            column_split = my_split, #Human_merge_metadata$Dataset_NEPC,
                            column_gap = unit(0.9,"in"),
                            column_title_gp = gpar(fontsize = ifelse(protein_set_name == "cis", 25, 60),
                                                   fontface = "bold",
                                                   col = c("#B71234", "#B71234", "#003884", "#003884", "#E09F3E", "#E09F3E")),
                            border_gp = gpar(col = "lightgray", lty = 1, lwd = 3, alpha = 0.6),
                            heatmap_legend_param = list(
                              title_gp = gpar(fontsize = ifelse(protein_set_name == "cis", 20, 35), fontface = "bold"),
                              labels_gp = gpar(fontsize = ifelse(protein_set_name == "cis", 20, 35))
                            )
    )
    ht_list_simpleAnnot = Human_merge_df_annot_cols_simpleAnnot %v% main_hm_gExpr %v% main_hm_pAct
    my_height = 24
    my_width = 70
  } else {
    ht_list_simpleAnnot = Human_merge_df_annot_cols_simpleAnnot %v% main_hm_pAct
    if(protein_set_name == "cis"){
      my_height = 12
      my_width = 70
    } else {
      my_height = 84
      my_width = 140
    }
  }
  
  
  
  p <- draw(ht_list_simpleAnnot,
            column_km = 1,
            merge_legend = TRUE,
            heatmap_legend_side = "bottom", 
            annotation_legend_side = "bottom")
  pdf(out_file, width = my_width , height = my_height)
  print(p)
  dev.off()
}

dir.create("experiments/heatmaps/")
# SAVE THE VIPER NES-Scores (From GExpr Z-Scores, dataset reference)
saveVIPERTopCISGenesHeatmap(Su2c_EC_pAct,
                            Su2c_WC_pAct,
                            Beltran_pAct,
                            Human_merge_metadata,
                            protein_set_name = "cis",
                            gExpr = FALSE,
                            out_file = "experiments/heatmaps/top_CIS_hits_VIPER_NES_datasetRef.pdf",
                            Su2c_EC_gExpr,
                            Su2c_WC_gExpr,
                            Beltran_gExpr)
saveVIPERTopCISGenesHeatmap(Su2c_EC_pAct,
                            Su2c_WC_pAct,
                            Beltran_pAct,
                            Human_merge_metadata,
                            protein_set_name = "mrs",
                            gExpr = FALSE,
                            out_file = "experiments/heatmaps/top_bottom_mrs_VIPER_NES_datasetRef.pdf")

# # SAVE THE GEXPR
# saveVIPERTopCISGenesHeatmap(Su2c_EC_gExpr,
#                             Su2c_WC_gExpr,
#                             Beltran_gExpr,
#                             Human_merge_metadata,
#                             protein_set_name = "cis",
#                             gExpr = TRUE,
#                             out_file = "experiments/heatmaps/top_CIS_hits_gExpr_datasetRef.pdf")
saveVIPERTopCISGenesHeatmap(Su2c_EC_gExpr,
                            Su2c_WC_gExpr,
                            Beltran_gExpr,
                            Human_merge_metadata,
                            protein_set_name = "mrs",
                            gExpr = TRUE,
                            out_file = "experiments/heatmaps/top_bottom_mrs_gExpr_datasetRef.pdf")

# SAVE THE VIPER Z-Scores (From GExpr Z-Scores, dataset reference) 
Su2c_EC_pAct_Z <- ZScoreTransform(Su2c_EC_pAct)
Su2c_WC_pAct_Z <- ZScoreTransform(Su2c_WC_pAct)
Beltran_pAct_Z <- ZScoreTransform(Beltran_pAct)
saveVIPERTopCISGenesHeatmap(Su2c_EC_pAct_Z,
                            Su2c_WC_pAct_Z,
                            Beltran_pAct_Z,
                            Human_merge_metadata,
                            protein_set_name = "cis",
                            gExpr = FALSE,
                            out_file = "experiments/heatmaps/top_CIS_hits_VIPER_ZScore_datasetRef.pdf",
                            Su2c_EC_gExpr,
                            Su2c_WC_gExpr,
                            Beltran_gExpr)
saveVIPERTopCISGenesHeatmap(Su2c_EC_pAct_Z,
                            Su2c_WC_pAct_Z,
                            Beltran_pAct_Z,
                            Human_merge_metadata,
                            protein_set_name = "mrs",
                            gExpr = FALSE,
                            out_file = "experiments/heatmaps/top_bottom_mrs_VIPER_ZScore_datasetRef.pdf")

# SAVE THE GExpr Z-Scores, global reference)
Human_gExpr_globalRef <- Su2c_EC_gExpr %>%
  mergeMatrices(Su2c_WC_gExpr) %>%
  mergeMatrices(Beltran_gExpr) %>%
  ZScoreTransform()
Su2c_EC_gExpr_globalRef <- Human_gExpr_globalRef[, colnames(Su2c_EC_gExpr)]
Su2c_WC_gExpr_globalRef <- Human_gExpr_globalRef[, colnames(Su2c_WC_gExpr)]
Beltran_gExpr_globalRef <- Human_gExpr_globalRef[, colnames(Beltran_gExpr)]
saveVIPERTopCISGenesHeatmap(Su2c_EC_gExpr_globalRef,
                            Su2c_WC_gExpr_globalRef,
                            Beltran_gExpr_globalRef,
                            Human_merge_metadata,
                            protein_set_name = "cis",
                            gExpr = TRUE,
                            out_file = "experiments/heatmaps/top_CIS_hits_gExpr_globalRef.pdf")
saveVIPERTopCISGenesHeatmap(Su2c_EC_gExpr_globalRef,
                            Su2c_WC_gExpr_globalRef,
                            Beltran_gExpr_globalRef,
                            Human_merge_metadata,
                            protein_set_name = "mrs",
                            gExpr = TRUE,
                            out_file = "experiments/heatmaps/top_bottom_mrs_gExpr_globalRef.pdf")

# SAVE THE VIPER NES-Scores (From GExpr Z-Scores, global reference)
su2cEC_regulon <- readRDS("data/SU2C_EastCoast/SU2CEastCoast_net_pruned.rds")
su2cWC_regulon <- readRDS("data/SU2C_WestCoast/SU2CWestCoast_net_pruned.rds")
Human_pAct_globalRef <- viper(Human_gExpr_globalRef,
                              list(su2cEC_regulon, su2cWC_regulon),
                              method = 'none',
                              eset.filter = FALSE,
                              mvws = 10)
Su2c_EC_pAct_globalRef <- Human_pAct_globalRef[, colnames(Su2c_EC_pAct)]
Su2c_WC_pAct_globalRef <- Human_pAct_globalRef[, colnames(Su2c_WC_pAct)]
Beltran_pAct_globalRef <- Human_pAct_globalRef[, colnames(Beltran_pAct)]
saveVIPERTopCISGenesHeatmap(Su2c_EC_pAct_globalRef,
                            Su2c_WC_pAct_globalRef,
                            Beltran_pAct_globalRef,
                            Human_merge_metadata,
                            protein_set_name = "cis",
                            gExpr = FALSE,
                            out_file = "experiments/heatmaps/top_CIS_hits_VIPER_NES_globalRef.pdf",
                            Su2c_EC_gExpr_globalRef,
                            Su2c_WC_gExpr_globalRef,
                            Beltran_gExpr_globalRef)
saveVIPERTopCISGenesHeatmap(Su2c_EC_pAct_globalRef,
                            Su2c_WC_pAct_globalRef,
                            Beltran_pAct_globalRef,
                            Human_merge_metadata,
                            protein_set_name = "mrs",
                            gExpr = FALSE,
                            out_file = "experiments/heatmaps/top_bottom_mrs_VIPER_NES_globalRef.pdf")

# SAVE THE VIPER Z-Scores (From GExpr Z-Scores, global reference)
Su2c_EC_pAct_globalRef_Z <- ZScoreTransform(Su2c_EC_pAct_globalRef)
Su2c_WC_pAct_globalRef_Z <- ZScoreTransform(Su2c_WC_pAct_globalRef)
Beltran_pAct_globalRef_Z <- ZScoreTransform(Beltran_pAct_globalRef)
saveVIPERTopCISGenesHeatmap(Su2c_EC_pAct_globalRef_Z,
                            Su2c_WC_pAct_globalRef_Z,
                            Beltran_pAct_globalRef_Z,
                            Human_merge_metadata,
                            protein_set_name = "cis",
                            gExpr = FALSE,
                            out_file = "experiments/heatmaps/top_CIS_hits_VIPER_ZScore_globalRef.pdf",
                            Su2c_EC_gExpr_globalRef,
                            Su2c_WC_gExpr_globalRef,
                            Beltran_gExpr_globalRef)
saveVIPERTopCISGenesHeatmap(Su2c_EC_pAct_globalRef_Z,
                            Su2c_WC_pAct_globalRef_Z,
                            Beltran_pAct_globalRef_Z,
                            Human_merge_metadata,
                            protein_set_name = "mrs",
                            gExpr = FALSE,
                            out_file = "experiments/heatmaps/top_bottom_mrs_VIPER_ZScore_globalRef.pdf")
