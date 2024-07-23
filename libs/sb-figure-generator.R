##
# SB Figures Generator
# ------------------------------------------------------------
# source("sources/sleeping-beauty/for-paper/sb-figure-generator.R")

# Libraries
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(viper))
suppressPackageStartupMessages(library(dendextend))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library(RColorBrewer))

# Functions
source("libs/tools/utils.R")
source("libs/tools/cross_species_utils.R")
source("libs/tools/net-enrichment.R")
source("libs/tools/atools-gsea.R")

#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################

getColumnsDendrogram_k2 <- function(mat){
  my_hclust <- hclust(d = as.dist(1-cor(mat)),method = "ward.D2")
  columns_dend = as.dendrogram(my_hclust)
  columns_dend = color_branches(columns_dend, k = 2) 		
  return(columns_dend)
}
getDFAnnotColsPathwaysHM <- function(df_annot.tibble, dna_rep_color_fun){
  df_annot.tibble <- df_annot.tibble %>% dplyr::select(genotype,
                                                       sb_status,
                                                       phenotype,
                                                       beltran_enrichment)
  df_annot.df <- as.data.frame(df_annot.tibble)
  df_annot_cols = HeatmapAnnotation( df = df_annot.df ,
                                     col = list( phenotype = c("NEPC"="steelblue","Adeno"="tan1") ,
                                                 genotype = c("NP"="khaki","NPp53"="orchid") ,
                                                 sb_status = c("activated"="red3","not-active"="gray75") ,
                                                 beltran_enrichment = dna_rep_color_fun )
  )
  return(df_annot_cols)
}

# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
# -------------------------------------------------------------------------------------------
# -------------------------------- GET MAIN HEATMAP FUNCTIONS  ------------------------------
# -------------------------------------------------------------------------------------------
# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&

# ------------------------------------- HELPER FUNCTIONS ------------------------------------

slice_mat_top_rows <- function(x,n=20) { print(nrow(x)) ; x <- x[ order( rowVars( x ) , decreasing = T )[1:n] , ] ; print(nrow(x)) ; return(x) }
getMatPathways <- function(mat, n_top_rows=20, PATHWAY = "HALLMARK"){
  mat_pathways <- mat[ grepl(PATHWAY , rownames(mat) ) , ]
  rownames(mat_pathways) <- gsub(paste0(PATHWAY, "_"),"",rownames(mat_pathways))
  mat_pathways <- slice_mat_top_rows(mat_pathways, n = n_top_rows)
  return(mat_pathways)
}
# getMatHallmarks <- function(mat, n_top_rows=20){
#   mat_hallmarks <- mat[ grepl("HALLMARK" , rownames(mat) ) , ]
#   rownames(mat_hallmarks) <- gsub("HALLMARK_","",rownames(mat_hallmarks))
#   mat_hallmarks <- slice_mat_top_rows(mat_hallmarks, n = n_top_rows)
#   return(mat_hallmarks)
# }
# getMatKegg <- function(mat, n_top_rows=20){
#   mat_kegg <- mat[ grepl("KEGG" , rownames(mat) ) , ]
#   rownames(mat_kegg) <- gsub("KEGG_","",rownames(mat_kegg))
#   mat_kegg <- slice_mat_top_rows(mat_kegg, n = n_top_rows)
#   return(mat_kegg)
# }
# getMatReactome <- function(mat, n_top_rows=20){
#   mat_reactome <- mat[ grepl("REACTOME" , rownames(mat) ) , ]
#   rownames(mat_reactome) <- gsub("REACTOME_","",rownames(mat_reactome))
#   mat_reactome <- slice_mat_top_rows(mat_reactome, n = n_top_rows)
#   return(mat_reactome)
# }
constructHallmarksHeatmap <- function(mat_hallmarks, col_fun, columns_dend){
  hallmarks_hm <- Heatmap( mat_hallmarks,
                           col = col_fun,
                           heatmap_legend_param = list(color_bar = "continuous"),
                           show_heatmap_legend = FALSE, 
                           cluster_columns = columns_dend,
                           name = "Cancer Hallmarks Pathway Enrichment",
                           row_names_gp = gpar(fontsize = 5),
                           column_names_gp = gpar(fontsize = 5),
                           row_title = "HALLMARKS", 
                           row_title_side = "left",
                           row_title_gp = gpar(fontsize = 10, fontface = "bold")
  )
  return(hallmarks_hm)
}

constructKeggHeatmap <- function(mat_kegg, col_fun, columns_dend){
  kegg_hm <- Heatmap( mat_kegg,
                      col = col_fun,
                      heatmap_legend_param = list(color_bar = "continuous"),
                      show_heatmap_legend = FALSE,
                      cluster_columns = columns_dend,
                      name = "KEGG Pathway Enrichment",
                      row_names_gp = gpar(fontsize = 5),
                      column_names_gp = gpar(fontsize = 5),
                      row_title = "KEGG", 
                      row_title_side = "left",
                      row_title_gp = gpar(fontsize = 10, fontface = "bold")											
  )
  return(kegg_hm)
}

constructReactomeHeatmap <- function(mat_reactome, col_fun, columns_dend){
  reactome_hm <- Heatmap( mat_reactome,
                          col = col_fun,
                          heatmap_legend_param = list(color_bar = "continuous"),
                          show_heatmap_legend = FALSE,
                          clustering_method_columns = "ward.D2",
                          name = "REACTOME Pathway Enrichment",
                          row_names_gp = gpar(fontsize = 5),
                          column_names_gp = gpar(fontsize = 5),
                          row_title = "REACTOME", 
                          row_title_side = "left",
                          row_title_gp = gpar(fontsize = 10, fontface = "bold")
  )
  return(reactome_hm)
}

# -------------------------------------- MAIN FUNCTIONS -------------------------------------
getHallmarksHeatmap <- function(mat, col_fun, columns_dend){
  mat_hallmarks <- getMatPathways(mat, PATHWAY = "HALLMARK")
  hallmarks_hm <- constructHallmarksHeatmap(mat_hallmarks, col_fun, columns_dend)
  return(hallmarks_hm)
}
getKeggHeatmap <- function(mat, col_fun, columns_dend){
  mat_kegg <- getMatPathways(mat, PATHWAY = "KEGG")
  kegg_hm <- constructKeggHeatmap(mat_kegg, col_fun, columns_dend)
  return(kegg_hm)
}
getReactomeHeatmap <- function(mat, col_fun, columns_dend){
  mat_reactome <- getMatPathways(mat, PATHWAY = "REACTOME")
  reactome_hm <- constructReactomeHeatmap(mat_reactome, col_fun, columns_dend)
  return(reactome_hm)
}

# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
# -------------------------------------------------------------------------------------------
# ---------------------------- GET PATHWAYS HEATMAP LIST FUNCTION ---------------------------
# -------------------------------------------------------------------------------------------
# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
getPathwaysHeatmapsList <- function(sb_data, filter_NP=TRUE){
  mat <- sb_data$pathways
  
  if(filter_NP == TRUE){
    mat <- mat[, sb_data$metadata$genotype == "NPp53"]
  }
  
  columns_dend <- getColumnsDendrogram_k2(mat)
  col_fun = colorRamp2( seq(-10 , 10 , length.out = 10 ) , RColorBrewer::brewer.pal(10,"Spectral") %>% rev() )
  
  hallmarks_hm <- getHallmarksHeatmap(mat, col_fun, columns_dend)
  kegg_hm <- getKeggHeatmap(mat, col_fun, columns_dend)
  reactome_hm <- getReactomeHeatmap(mat, col_fun, columns_dend)
  pws_hms_list <- list(hallmarks_hm, kegg_hm, reactome_hm)
  names(pws_hms_list) <- c("hallmarks_hm", "kegg_hm", "reactome_hm")
  return(pws_hms_list)
}


# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
# -------------------------------------------------------------------------------------------
# ----------------------------- SAVE PATHWAYS HEATMAP FUNCTIONS -----------------------------
# -------------------------------------------------------------------------------------------
# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&

# ------------------------------------- HELPER FUNCTIONS ------------------------------------

savePathwaysHeatmap <- function(df_annot_cols,
                                hallmarks_hm,
                                kegg_hm,
                                reactome_hm,
                                dna_rep_color_fun){
  ht_list = df_annot_cols %v% hallmarks_hm %v% kegg_hm %v% reactome_hm
  lgd_pathway_enrichment = getPathwayEnrichmentLegend()
  lgd_beltran_enrichment = getBeltranEnrichmentLegend(dna_rep_color_fun)
  pd = packLegend(lgd_beltran_enrichment,
                  lgd_pathway_enrichment,
                  direction = "hor",
                  max_width = unit(4, "in") )
  p <- draw(ht_list,
            column_km = 1, 
            annotation_legend_list = pd,
            heatmap_legend_side = "bottom", 
            annotation_legend_side = "bottom" )
  
  pdf( file.path( "experiments/sb-figure-generator/reports/sb-pathways-heatmap.pdf") , width = 6 , height = 8 )
  print(p)
  dev.off()
}

# ---------------------------------- MAIN FUNCTION ----------------------------------
savePathwaysHeatmapsAsPDFs <- function(sb_data, pws_hms_list){
  df_annot.tibble <- sb_data$metadata
  mat <- sb_data$pathways
  stopifnot( identical( colnames(mat) , df_annot.tibble$sample_id ) )
  
  dna_rep_color_fun = colorRamp2(c(-5, 0, 5), c("blue", "white", "red"))
  
  df_annot_cols <- getDFAnnotColsPathwaysHM(df_annot.tibble, dna_rep_color_fun)
  
  savePathwaysHeatmap(df_annot_cols,
                      pws_hms_list[["hallmarks_hm"]],
                      pws_hms_list[["kegg_hm"]],
                      pws_hms_list[["reactome_hm"]],
                      dna_rep_color_fun)
}


#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################

#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################

# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
# -----------------------------------------------------------------------------------
# -------------------------------- LOAD DATA FUNCTIONS ------------------------------
# -----------------------------------------------------------------------------------
# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
loadCISTablePerLibrary <- function(){
  my_filename <- file.path("experiments/integrative-analysis/processed_data/cis-table-per-library.rds")
  cis_table.per_library <- readRDS(my_filename)
  return(cis_table.per_library)
}
loadRNASeqTable <- function(){
  my_filename <- file.path("experiments/all-samples-sb-rna-seq-analysis/reports/sb-rna-seq-samples-table.xlsx")
  rna_seq_table <- readxl::read_excel(my_filename,sheet = 1)
  colnames(rna_seq_table) <- make.names(colnames(rna_seq_table))
  
  index <- !is.na(rna_seq_table$sample_id)
  rna_seq_table <- rna_seq_table[index,]
  return(rna_seq_table)
}
loadDNASeqTable <- function(){
  my_filename <- file.path("data/gene-expression-from-sb-mouse-models/RNAseq and CIS list Sep2020.xlsx")
  dna_seq_table <- readxl::read_excel(my_filename,sheet = 3)
  colnames(dna_seq_table) <- make.names(colnames(dna_seq_table))
  
  dna_seq_table$Tumor.ID_mod <- gsub( "(.*) (.*)" , "\\1" , dna_seq_table$Tumor.ID )
  return(dna_seq_table)
}

# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
# -----------------------------------------------------------------------------------
# ------------------------------ COMBINE DATA FUNCTIONS -----------------------------
# -----------------------------------------------------------------------------------
# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
joinRNA_DNATables <- function(rna_seq_table , dna_seq_table){
  joined_tables <- left_join( rna_seq_table , dna_seq_table , by = c("mouse.ID"="Tumor.ID_mod"))
  stopifnot(identical( ifelse( joined_tables$has_dna_paired_samples == "Yes" , TRUE , FALSE ) , !is.na(joined_tables$Tumor.ID) ))
  
  index <- !is.na(joined_tables$DNA.Sample.ID)
  joined_tables <- joined_tables[index,]
  
  colnames(joined_tables)
  joined_tables <- joined_tables %>% dplyr::select(sample_id,mouse.ID,DNA.Sample.ID)	
  return(joined_tables)
}
joinCIS_RNA_DNATables <- function(cis_table.per_library, rna_seq_table , dna_seq_table){
  joined_tables <- joinRNA_DNATables(rna_seq_table , dna_seq_table)
  cis_table.per_library <- left_join( cis_table.per_library , joined_tables , by = c("library_name"="DNA.Sample.ID") )
  return(cis_table.per_library)
}
getRNASeqSamplesThatHaveDNASeqSamples <- function(cis_table.per_library, rna_seq_table){
  joined_tables <- joinRNA_DNATables(rna_seq_table, dna_seq_table)
  rnaseq_samples_that_have_dnaseq_samples <- unique(joined_tables$sample_id)
  # length(rnaseq_samples_that_have_dnaseq_samples)
  rnaseq_samples_that_have_dnaseq_samples.table <- table(joined_tables$sample_id)
  return(rnaseq_samples_that_have_dnaseq_samples.table)
}
filterCISTableToTopHits <- function(cis_table.per_library, my_top_hits_table){
  cis_table.per_library <- cis_table.per_library %>% 
    filter( mouse_gene_name %in% my_top_hits_table$gene ) %>%
    filter( !is.na(sample_id) )
  return(cis_table.per_library)
}

# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
# -----------------------------------------------------------------------------------
# ------------------------- CONSTRUCT GRAPH MATRIX FUNCTIONS ------------------------
# -----------------------------------------------------------------------------------
# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&

# --------------------------------- HELPER FUNCTION ---------------------------------
constructInitialGraphMatrix <- function(cis_table.per_library, sb_data){
  my_metadata <- sb_data$metadata
  
  tmp <- inner_join( my_metadata , cis_table.per_library , by = c("sample_id"="sample_id") )
  tmp <- tmp %>% group_by(sample_id) %>% summarise(genes = paste(mouse_gene_name, collapse=","))
  
  tmp <- tmp %>% separate_rows(genes,sep = ",")
  
  
  g <- graph_from_data_frame(tmp %>% as.data.frame())
  V(g)$type <- grepl("CMZ*" , names(V(g)) )
  g_mat <- as_incidence_matrix( g )
  
  return(g_mat)
}
addRemainingSamplesToGraphMatrix <- function(g_mat, sb_data){
  df_annot.tibble <- sb_data$metadata
  
  remaining_samples <- setdiff( df_annot.tibble$sample_id , colnames(g_mat) )
  
  r_mat <- matrix( 0L,
                   nrow = nrow(g_mat),
                   ncol = length(remaining_samples),
                   dimnames = list(rownames(g_mat),remaining_samples) )
  
  g_mat <- cbind(g_mat, r_mat)
  
  index <- match( df_annot.tibble$sample_id, colnames(g_mat) )
  g_mat <- g_mat[,index]
  stopifnot(identical(colnames(g_mat),df_annot.tibble$sample_id))
  
  return(g_mat)
}

# ---------------------------------- MAIN FUNCTION ----------------------------------
constructCompleteGraphMatrix <- function(cis_table.per_library, sb_data){
  g_mat <- cis_table.per_library %>%
    constructInitialGraphMatrix(sb_data) %>%
    addRemainingSamplesToGraphMatrix(sb_data)
  return(g_mat)
}

# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
# -----------------------------------------------------------------------------------
# ------------------------- CONSTRUCT GRAPH MATRIX MAIN FUNC ------------------------
# -----------------------------------------------------------------------------------
# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&

combineDataToConstructGraphMatrixWithCISHits <- function(cis_table.per_library,
                                                         rna_seq_table,
                                                         dna_seq_table,
                                                         my_top_hits_table,
                                                         sb_data) {
  # Join Data
  cis_table.per_library <- joinCIS_RNA_DNATables(cis_table.per_library, rna_seq_table , dna_seq_table)
  my_filename <- file.path("experiments/sb-figure-generator/processed_data/cis-table-per-library-with-dna-sample-id.rds")
  saveRDS(cis_table.per_library, my_filename)
  # View( cis_table.per_library %>% 
  #         group_by(mouse_gene_name,phenotype) %>% tally() %>%
  #         arrange(desc(n),mouse_gene_name,phenotype) %>%
  #         group_by(mouse_gene_name) %>% summarize(mouse_gene_name,phenotype,fold=n/n[phenotype=="adeno"])
  # )
  
  # ---
  cis_table.per_library <- filterCISTableToTopHits(cis_table.per_library, my_top_hits_table)
  g_mat <- constructCompleteGraphMatrix(cis_table.per_library, sb_data)
  return(g_mat)
}

#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################

#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################

getSigList <- function(path_to_sig_csv){
  read_csv(path_to_sig_csv, col_names = FALSE)$X1
}
getMatWithOnlyTFAndCoTF <- function(mat, sig_list){
  index <- !(rownames(mat) %in% sig_list)
  mat <- mat[index,]
  return(mat)
}
getTopNProteinsFromStoufferSignature <- function(mat, n_top = 15){
  mat.stouffered <- doStoufferFromClusters( mat , factor(kmeans(t(mat), centers = 3,nstart = 15)$cluster) )
  index <- names(sort(mat.stouffered["Ar",],decreasing = T))
  mat.stouffered <- mat.stouffered[ , index ]
  selected_proteins <- Reduce( union , apply( mat.stouffered , 2 , function(x) { x <- sort(x,decreasing = T) ; names(c(head(x,n_top))) } ) )
  return(selected_proteins)
}
getColumnsDendrogram_k3 <- function(mat){
  my_hclust <- hclust(d = as.dist(1-cor(mat,method = "spe")),method = "ward.D2")
  columns_dend = as.dendrogram(my_hclust)
  columns_dend = color_branches(columns_dend, k = 3)
  return(columns_dend)
}

# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
# -------------------------------------------------------------------------------------------
# -------------------- COMPLEXHEATMAP: GET DF ANNOTATION COLUMN FUNCTION --------------------
# -------------------------------------------------------------------------------------------
# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&

# ------------------------------------- HELPER FUNCTIONS ------------------------------------
getDfAnnotDF <- function(my_metadata, mat, beltran_only = FALSE){
  if(beltran_only){
    df_annot.tibble <- my_metadata %>% dplyr::select(rna_dna_matched_samples,beltran_enrichment)
  } else {
    df_annot.tibble <- my_metadata %>% dplyr::select(phenotype,genotype,sb_status,qPCR.Syn,qPCR.cga)
  }
  df_annot.df <- as.data.frame(df_annot.tibble)
  rownames(df_annot.df) <- colnames(mat)
  return(df_annot.df)
}
# -------------------------------------- MAIN FUNCTIONS -------------------------------------
getDfAnnotCols <- function(my_metadata, mat, qPCR_color_fun){
  df_annot.df <- getDfAnnotDF(my_metadata, mat)
  df_annot_cols = HeatmapAnnotation( df = df_annot.df ,
                                     col = list( phenotype = c("NEPC"="steelblue","Adeno"="tan1") ,
                                                 genotype = c("NP"="khaki","NPp53"="orchid") ,
                                                 sb_status = c("activated"="red3","not-active"="gray75") ,
                                                 qPCR.Syn = qPCR_color_fun ,
                                                 qPCR.cga = qPCR_color_fun 
                                     ) , show_legend = TRUE , na_col = "gray95"
  ) 
  return(df_annot_cols)
}
getDfAnnotColsBeltranOnly <- function(my_metadata, mat, dna_rep_color_fun){
  df_annot.df_beltran_only <- getDfAnnotDF(my_metadata, mat, beltran_only = TRUE)
  df_annot_cols_beltran_only = columnAnnotation( df = df_annot.df_beltran_only[,"beltran_enrichment",drop=FALSE] ,
                                                 col = list(
                                                   beltran_enrichment = dna_rep_color_fun
                                                 ) , show_legend = TRUE , na_col = "gray95"
  )
  return(df_annot_cols_beltran_only)
}

# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
# -----------------------------------------------------------------------------------
# ----------------------------- GET MAIN HEATMAPS FUNCS -----------------------------
# -----------------------------------------------------------------------------------
# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
getCISHeatmap <- function(g_mat, columns_dend){
  cis_hm <- Heatmap(
    g_mat ,
    cell_fun = function(j, i, x, y, w, h, col) {
      # add text to each grid
      grid.text(round(g_mat[i, j], 1), x, y, gp = gpar(fontsize = 5, col = "gray90"))
    },
    col = circlize::colorRamp2(c(0, 1, 3), c("gray98", "gray5", "gray1")),
    heatmap_legend_param = list(color_bar = "continuous"),
    show_heatmap_legend = FALSE,
    cluster_columns = columns_dend,
    name = "CIS-associated Genes",
    rect_gp = gpar(col = "white", lwd = 1),
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 8),
    row_title = "CIS-associated Gene",
    row_title_side = "left",
    row_title_gp = gpar(fontsize = 10, fontface = "bold")
  )	
  return(cis_hm)
}
getNESHeatmap <- function(mat, columns_dend){
  nes_hm <- Heatmap( mat,
                     col = colorRamp2(c(-3, 0, 3), c("deepskyblue3", "white", "brown3")),
                     heatmap_legend_param = list(color_bar = "continuous"),
                     show_heatmap_legend = FALSE, 
                     cluster_rows = FALSE,
                     cluster_columns = columns_dend,
                     name = "TF Activity",
                     row_names_gp = gpar(fontsize = 10),
                     column_names_gp = gpar(fontsize = 10),
                     row_title = "TF Activity", 
                     row_title_side = "left",
                     row_title_gp = gpar(fontsize = 12, fontface = "bold"),
                     use_raster = FALSE
  )
  return(nes_hm)
}
getNESHeatmapWithKnownMarkers <- function(mat_known, columns_dend){
  nes_known_hm <- Heatmap( mat_known,
                           col = colorRamp2(c(-7, 0, 7), c("deepskyblue3", "white", "brown3")),
                           heatmap_legend_param = list(color_bar = "continuous"),
                           show_heatmap_legend = FALSE, 
                           cluster_rows = FALSE,
                           cluster_columns = columns_dend,
                           name = "TF Activity",  
                           row_names_gp = gpar(fontsize = 10),
                           column_names_gp = gpar(fontsize = 10),
                           row_title = "TF Activity", 
                           row_title_side = "left",
                           row_title_gp = gpar(fontsize = 12, fontface = "bold"),
                           use_raster = FALSE
  )
  return(nes_known_hm)
}

# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
# -----------------------------------------------------------------------------------
# -------------------------------- SAVE HEATMAP FUNCS -------------------------------
# -----------------------------------------------------------------------------------
# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&

# -------------------------------- HELPER FUNCTIONS ---------------------------------
getProteinActivityScoreEnrichmentLegend <- function(col_fun){
  lgd_pas_enrichment = Legend( col_fun = col_fun,
                               title = "Protein Activity Score",
                               at = c(-10 , -1.645 , 0, +1.645 , +10 ),
                               labels = c("-10","Lower Threshold", "Zero", "Higher Threshold","+10"),
                               direction = "horizontal",
                               title_position = "topcenter")	
  return(lgd_pas_enrichment)
}
getBeltranEnrichmentLegend <- function(dna_rep_color_fun){
  lgd_beltran_enrichment = Legend( col_fun = dna_rep_color_fun,
                                   title = "Beltran Enrichment Score",
                                   at = c(-10 , -1.645 , 0, +1.645 , +10 ) ,
                                   labels = c("-10","Lower Threshold", "Zero", "Higher Threshold","+10"),
                                   direction = "horizontal",
                                   title_position = "topcenter")
  return(lgd_beltran_enrichment)
}
getPathwayEnrichmentLegend <- function(){
  col_fun = colorRamp2( seq(-10 , 10 , length.out = 10 ), RColorBrewer::brewer.pal(10,"Spectral") %>% rev() )
  lgd_pathway_enrichment = Legend( col_fun = col_fun,
                                   title = "Pathway Enrichment Score",
                                   at = c(-10 , -1.645 , 0, +1.645 , +10 ),
                                   labels = c("-10","Lower Threshold", "Zero", "Higher Threshold","+10"),
                                   direction = "horizontal",
                                   title_position = "topcenter")	
  return(lgd_pathway_enrichment)
}
  


# --------------------------------- MAIN FUNCTIONS ----------------------------------
saveVIPERPlusCISHeatmap <- function(df_annot_cols,
                                    cis_hm,
                                    df_annot_cols_beltran_only,
                                    nes_hm,
                                    dna_rep_color_fun,
                                    col_fun){
  ht_list = df_annot_cols %v% cis_hm %v% df_annot_cols_beltran_only %v% nes_hm
  lgd_beltran_enrichment <- getBeltranEnrichmentLegend(dna_rep_color_fun)
  lgd_pas_enrichment <- getProteinActivityScoreEnrichmentLegend(col_fun)
  
  pd = packLegend(lgd_beltran_enrichment, lgd_pas_enrichment,direction = "hor",max_width = unit(4, "in") )
  p <- draw(ht_list, column_km = 1 , 
            annotation_legend_list = pd ,
            heatmap_legend_side = "bottom", 
            annotation_legend_side = "bottom" )
  
  pdf( file.path( "experiments/sb-figure-generator/reports/sb-cis-heatmap.pdf") , width = 8 , height = 12 )
  print(p)
  dev.off()
}
saveTFsActivityPathwaysHeatmap <- function(df_annot_cols,
                                           cis_hm,
                                           df_annot_cols_beltran_only,
                                           nes_hm,
                                           hallmarks_hm,
                                           kegg_hm,
                                           reactome_hm,
                                           dna_rep_color_fun,
                                           col_fun){
  ht_list = df_annot_cols %v% cis_hm %v% df_annot_cols_beltran_only %v% nes_hm %v% hallmarks_hm %v% kegg_hm %v% reactome_hm
  
  lgd_beltran_enrichment <- getBeltranEnrichmentLegend(dna_rep_color_fun)
  lgd_pas_enrichment <- getProteinActivityScoreEnrichmentLegend(col_fun)
  lgd_pathway_enrichment <- getPathwayEnrichmentLegend()
  
  pd = packLegend(lgd_beltran_enrichment,
                  lgd_pathway_enrichment,
                  lgd_pas_enrichment,direction = "hor",
                  max_width = unit(4, "in") )
  
  p <- draw(ht_list, column_km = 1 , 
            annotation_legend_list = pd ,
            heatmap_legend_side = "bottom", 
            annotation_legend_side = "bottom" )
  pdf( file.path( "experiments/sb-figure-generator/reports/sb-cis-nes-pathways-heatmap.pdf") , width = 8 , height = 18 )
  print(p)
  dev.off()  			
}
saveKnownTFsActivityPathwaysHeatmap<- function(df_annot_cols,
                                               df_annot_cols_beltran_only,
                                               nes_known_hm,
                                               hallmarks_hm,
                                               dna_rep_color_fun,
                                               col_fun){
  ht_list = df_annot_cols %v% df_annot_cols_beltran_only %v% nes_known_hm %v% hallmarks_hm
  
  lgd_beltran_enrichment <- getBeltranEnrichmentLegend(dna_rep_color_fun)
  lgd_pas_enrichment <- getProteinActivityScoreEnrichmentLegend(col_fun)
  lgd_pathway_enrichment <- getPathwayEnrichmentLegend()
  
  pd = packLegend(lgd_beltran_enrichment,
                  lgd_pathway_enrichment,
                  lgd_pas_enrichment,
                  direction = "hor",
                  max_width = unit(4, "in") )
  
  p <- draw(ht_list, column_km = 1 , 
            annotation_legend_list = pd ,
            heatmap_legend_side = "bottom", 
            annotation_legend_side = "bottom" )
  pdf( file.path( "experiments/sb-figure-generator/reports/sb-kwown-nes-pathways-heatmap.pdf") , width = 7 , height = 8 )
  print(p)
  dev.off()  			
}

# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
# -----------------------------------------------------------------------------------
# ----------------------- PRINT VIPER HEATMAP PLUS CIS FUNCTION ---------------------
# -----------------------------------------------------------------------------------
# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
printVIPERHeatmapPlusCIS <- function(sb_data,
                                     g_mat,
                                     sig_list,
                                     cis_table.per_library,
                                     rna_seq_table,
                                     hallmarks_hm,
                                     kegg_hm,
                                     reactome_hm,
                                     filter_NP=TRUE){
  my_metadata <- sb_data$metadata
  mat <- sb_data$vpmat
  
  print_msg_warn(">>> >> Selecting only TFs and coTFs ...")	
  mat <- getMatWithOnlyTFAndCoTF(mat, sig_list)
  selected_proteins <- getTopNProteinsFromStoufferSignature(mat)
  
  if(filter_NP == TRUE){
    my_metadata <- sb_data$metadata %>%
      filter(genotype == "NPp53")
    mat <- sb_data$vpmat[, sb_data$metadata$genotype == "NPp53"]
    g_mat <- g_mat[, sb_data$metadata$genotype == "NPp53"]
  }
  
  columns_dend <- getColumnsDendrogram_k3(mat)
  mat_known <- mat[ c("Ar","Ascl1","Foxa2","Foxm1","Cenpf") , ]
  mat <- mat[selected_proteins,]
  
  # Add RNA DNA Matched Samples To Metadata
  rnaseq_samples_that_have_dnaseq_samples.table <-
    getRNASeqSamplesThatHaveDNASeqSamples(cis_table.per_library, rna_seq_table)
  
  my_metadata$rna_dna_matched_samples <- rnaseq_samples_that_have_dnaseq_samples.table[ match(my_metadata$sample_id,names(rnaseq_samples_that_have_dnaseq_samples.table)) ]
  
  # Color Functions
  qPCR_color_fun =  colorRamp2( c(0.01,1,2,5,10,20,50,100,200,300) , viridis::viridis(10) )
  dna_rep_color_fun = colorRamp2( seq(-5 , 5 , length.out = 10 ) , RColorBrewer::brewer.pal(10,"Spectral") %>% rev() )
  col_fun = colorRamp2(c(-3, 0, 3), c("deepskyblue3", "white", "brown3"))
  
  # DF Annot Cols
  df_annot_cols <- getDfAnnotCols(my_metadata, mat, qPCR_color_fun)
  df_annot_cols_beltran_only <- getDfAnnotColsBeltranOnly(my_metadata, mat, dna_rep_color_fun)
  
  # Get Heatmaps
  cis_hm <- getCISHeatmap(g_mat, columns_dend)
  nes_hm <- getNESHeatmap(mat, columns_dend)
  nes_known_hm <- getNESHeatmapWithKnownMarkers(mat_known, columns_dend)
  
  # Save VIPER-CIS Heatmap
  saveVIPERPlusCISHeatmap(
    df_annot_cols,
    cis_hm,
    df_annot_cols_beltran_only,
    nes_hm,
    dna_rep_color_fun,
    col_fun
  )
  
  print_msg_info(">>> >> Printing TFs activity Pathways")
  saveTFsActivityPathwaysHeatmap(
    df_annot_cols,
    cis_hm,
    df_annot_cols_beltran_only,
    nes_hm,
    hallmarks_hm,
    kegg_hm,
    reactome_hm,
    dna_rep_color_fun,
    col_fun
  )
  
  print_msg_info(">>> >> Printing only Known TFs activity")
  saveKnownTFsActivityPathwaysHeatmap(
    df_annot_cols,
    df_annot_cols_beltran_only,
    nes_known_hm,
    hallmarks_hm,
    dna_rep_color_fun,
    col_fun
  )
}

#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################

# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
# -----------------------------------------------------------------------------------
# ---------------------------- SAVE CIS GEXPR HEATMAP FUNC --------------------------
# -----------------------------------------------------------------------------------
# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
saveCISGeneExpressionHeatmap <- function(sb_data, cis_table.per_library){
  cis_list <- cis_table.per_library$mouse_gene_name %>% unique()
  print(length(cis_list))
  
  my_metadata <- sb_data$metadata
  mat <- sb_data$raw_counts %>% edgeR::cpm(log=TRUE)
  
  df_annot.tibble <- my_metadata %>% dplyr::select(genotype,sb_status,phenotype)
  df_annot.df <- as.data.frame(df_annot.tibble)
  df_annot_cols = HeatmapAnnotation( df = df_annot.df ,
                                     col = list( phenotype = c("NEPC"="steelblue","Adeno"="tan1") ,
                                                 genotype = c("NP"="khaki","NPp53"="orchid") ,
                                                 sb_status = c("activated"="red3","not-active"="gray75")
                                     )
  )
  
  mat <- mat[ rownames(mat) %in% cis_list , ]
  mat <- rankNorm(mat)
  
  col_fun = colorRamp2( seq(-5 , 5 , length.out = 10 ) ,  brewer.pal(10,"PuOr") %>% rev() )
  gex_hm <- Heatmap( mat ,
                     col = col_fun ,
                     heatmap_legend_param = list(color_bar = "continuous") ,	
                     show_heatmap_legend = TRUE ,
                     cluster_rows = TRUE ,
                     cluster_columns = hclust(as.dist(1-cor(mat,method = "spe")),"ward.D2") ,
                     name = "Gene Expression (z-score)" ,
                     show_row_names = FALSE , show_column_names = FALSE ,
                     row_title = "Expression of CIS-associated Genes",
                     row_title_side = "left" , row_title_gp = gpar(fontsize = 10, fontface = "bold")
  )
  
  ht_list = df_annot_cols %v% gex_hm
  
  p <- draw(ht_list, column_km = 1 )
  pdf( file.path( "experiments/sb-figure-generator/reports/sb-cis-gene-expression-heatmap.pdf") ,
       width = 7 ,
       height = 5 )
  print(p)
  dev.off()
}

#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################

saveGSEAAnalysisOfBeltranAsPDF <- function(sb_data){
  tumor_sample <- sb_data$vpmat_from_edgeR_signatures[,"NPp53_NEPC_effect_on_SB_activated"]
  tumor_sample <- tumor_sample %>% mouse_to_human()
  gene_set <- beltran.net.regulon.obj[[1]]$tfmode # from sources/net-enrichment.R
  
  pdf(  file = file.path( "experiments/sb-figure-generator/reports/GSEA-SB-NEPC-over-Beltran-Gene-Set.pdf" ),
        width = 5,
        height = 4 )
  my_gsea <- gsea(
    signature = tumor_sample,
    geneset = gene_set,
    twoTails = TRUE,
    pout = TRUE,
    per = 1000,
    xlab = "Beltran, et al (2016) GS over SB-derived NEPC Signature",
    colSig = c(0.6, 1, 1, 0.85),
    lwd = 2
  )
  dev.off()	
}

#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################

loadBeltranGeneSet <- function(){
  filename <- "data/Beltran-2016/processed-data/novel-NEPC-gene-set-from-beltran-2016-data-interation-2.rds"
  beltran_new_regulon.obj <- readRDS(filename)
  gene_set <- beltran_new_regulon.obj[[1]]$tfmode 
  return(gene_set)
}
saveGSEAAnalysisOfNovelGeneSetAsPDF <- function(sb_data){
  tumor_sample <- sb_data$vpmat_from_edgeR_signatures[,"NPp53_NEPC_effect_on_SB_activated"]
  tumor_sample <- tumor_sample %>% mouse_to_human()
  gene_set <- loadBeltranGeneSet()
  
  pdf(  file = file.path( "experiments/sb-figure-generator/reports/GSEA-SB-NEPC-over-My-Novel-Beltran-Gene-Set.pdf" ),
        width = 5,
        height = 4 )
  my_gsea <- gsea(
    signature = tumor_sample,
    geneset = gene_set,
    twoTails = TRUE,
    pout = TRUE,
    per = 1000,
    xlab = "New Human NEPC Gene Set over SB-derived NEPC Signature",
    colSig = c(0.6, 1, 1, 0.85),
    lwd = 2
  )
  dev.off()	
}

#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################

#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################

# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
# -----------------------------------------------------------------------------------
# ---------------------- PRINT VIPER HEATMAP MINUS CIS FUNCTION ---------------------
# -----------------------------------------------------------------------------------
# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&

# --------------------------------- HELPER FUNCTION ---------------------------------
saveVIPERMinusCISHeatmap <- function(df_annot_cols,
                                     df_annot_cols_beltran_only,
                                     nes_hm,
                                     dna_rep_color_fun,
                                     col_fun,
                                     filter_NP,
                                     nes_known_hm) {
  ht_list = df_annot_cols %v% df_annot_cols_beltran_only %v% nes_hm %v% nes_known_hm
  lgd_beltran_enrichment <- getBeltranEnrichmentLegend(dna_rep_color_fun)
  lgd_pas_enrichment <- getProteinActivityScoreEnrichmentLegend(col_fun)
  
  pd = packLegend(lgd_beltran_enrichment, lgd_pas_enrichment,direction = "hor",max_width = unit(4, "in") )
  p <- draw(ht_list, column_km = 1 , 
            annotation_legend_list = pd ,
            heatmap_legend_side = "bottom", 
            annotation_legend_side = "bottom"#,
            #use_raster = FALSE
            )
  
  pdf( paste0( "experiments/sb-figure-generator/reports/sb-viper-",
               ifelse(filter_NP, "noNP-", ""),
               "heatmap.pdf") , width = 8 , height = 16 )
  print(p)
  dev.off()
}

# ---------------------------------- MAIN FUNCTION ----------------------------------
printVIPERHeatmapMinusCIS <- function(sb_data,
                                     g_mat,
                                     sig_list,
                                     cis_table.per_library,
                                     rna_seq_table,
                                     hallmarks_hm,
                                     kegg_hm,
                                     reactome_hm,
                                     filter_NP=TRUE){
  my_metadata <- sb_data$metadata
  mat <- sb_data$vpmat
  
  print_msg_warn(">>> >> Selecting only TFs and coTFs ...")	
  mat <- getMatWithOnlyTFAndCoTF(mat, sig_list)
  selected_proteins <- getTopNProteinsFromStoufferSignature(mat)
  
  columns_dend <- getColumnsDendrogram_k3(mat)
  nepc_genes_to_plot <- c("Ar",
                          "Foxa1",
                          "Nkx3-1",
                          "Foxa2",
                          "Nr3c1",
                          "Nsd2",
                          "Foxm1",
                          "Cenpf")
  mat_known <- mat[ nepc_genes_to_plot , ]
  mat <- mat[selected_proteins,]
  
  if(filter_NP == TRUE){
    my_metadata <- my_metadata %>%
      filter(genotype == "NPp53")
    mat <- mat[, my_metadata$sample_id] #[, sb_data$metadata$genotype == "NPp53"]
    mat_known <- mat_known[, my_metadata$sample_id]
    g_mat <- g_mat[, my_metadata$sample_id] #[, sb_data$metadata$genotype == "NPp53"]
    
    samples_to_remove <- sb_data$metadata$sample_id[sb_data$metadata$genotype=="NP"]
    columns_dend <- dendextend::prune(columns_dend, leaves = samples_to_remove)
  }
  
  # Add RNA DNA Matched Samples To Metadata
  rnaseq_samples_that_have_dnaseq_samples.table <- getRNASeqSamplesThatHaveDNASeqSamples(
    cis_table.per_library,
    rna_seq_table
  )
  
  my_metadata$rna_dna_matched_samples <- rnaseq_samples_that_have_dnaseq_samples.table[ match(my_metadata$sample_id, names(rnaseq_samples_that_have_dnaseq_samples.table)) ]
  
  # Color Functions
  qPCR_color_fun =  colorRamp2( c(0.01,1,2,5,10,20,50,100,200,300) , viridis::viridis(10) )
  dna_rep_color_fun = colorRamp2( seq(-5 , 5 , length.out = 10 ) ,
                                  RColorBrewer::brewer.pal(10,"Spectral") %>% rev() )
  col_fun = colorRamp2(c(-3, 0, 3), c("deepskyblue3", "white", "brown3"))
  
  # DF Annot Cols
  df_annot_cols <- getDfAnnotCols(my_metadata, mat, qPCR_color_fun)
  df_annot_cols_beltran_only <- getDfAnnotColsBeltranOnly(my_metadata, mat, dna_rep_color_fun)
  
  # Get Heatmaps
  nes_hm <- getNESHeatmap(mat, columns_dend)
  nes_known_hm <- getNESHeatmapWithKnownMarkers(mat_known, columns_dend)
  
  # Save VIPER-CIS Heatmap
  saveVIPERMinusCISHeatmap(
    df_annot_cols,
    df_annot_cols_beltran_only,
    nes_hm,
    dna_rep_color_fun,
    col_fun,
    filter_NP,
    nes_known_hm
  )
}

#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################

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

# Create Directories
dir.create("experiments/sb-figure-generator/", showWarnings = FALSE)
dir.create("experiments/sb-figure-generator/processed_data/", showWarnings = FALSE)
dir.create("experiments/sb-figure-generator/reports/", showWarnings = FALSE)

# Load Data
sb_data <- readRDS(file.path("experiments/all-samples-sb-rna-seq-analysis/processed_data/sb-rna-seq-data-analysis.rds"))

## Printing Pathway Heatmap ----
print_msg_info(">>> Printing Pathway Heatmap ...")
pws_hms_list <- getPathwaysHeatmapsList(sb_data, filter_NP=FALSE)
savePathwaysHeatmapsAsPDFs(sb_data, pws_hms_list)

## Generating Matrix with CIS hits for heatmap visualization ----
print_msg_info(">>> Generating Matrix with CIS hits for heatmap visualization ...")
cis_table.per_library <- loadCISTablePerLibrary()
rna_seq_table <- loadRNASeqTable()
dna_seq_table <- loadDNASeqTable()
my_top_hits_table <-
  readxl::read_excel("experiments/integrative-analysis/reports/nepc-signatures-table-all-3-conditions.xlsx")
g_mat <- combineDataToConstructGraphMatrixWithCISHits(cis_table.per_library,
                                                      rna_seq_table,
                                                      dna_seq_table,
                                                      my_top_hits_table,
                                                      sb_data)
## Loading Dataset ----
print_msg_info(">>> Loading Dataset ...")
# vpmat_vpsig <- readRDS("experiments/all-samples-sb-rna-seq-analysis/sb_vpmat_vpsig_npRef.rds")
sig_list <- getSigList("data/regulators-lists/sig.csv")

## Printing CIS hits on gene expression heatmap ----
print_msg_info(">>> >> Printing CIS hits on gene expression heatmap")
saveCISGeneExpressionHeatmap(sb_data, cis_table.per_library)

## Printing VIPER Heatmap + CIS ----
print_msg_info(">>> Printing VIPER Heatmap + CIS ...")
printVIPERHeatmapPlusCIS(sb_data,
                         g_mat,
                         sig_list,
                         cis_table.per_library,
                         rna_seq_table,
                         pws_hms_list[["hallmarks_hm"]],
                         pws_hms_list[["kegg_hm"]],
                         pws_hms_list[["reactome_hm"]],
                         filter_NP=FALSE)

## GSEA Analysis of Beltran, et al. Gene Set over SB NEPC signature ----
print_msg_info(">>> GSEA Analysis of Beltran, et al. Gene Set over SB NEPC signature")
saveGSEAAnalysisOfBeltranAsPDF(sb_data)

## GSEA Analysis of My Novel NEPC Gene Set identified using Beltran, et al., over SB NEPC signature ----
print_msg_info(">>> GSEA Analysis of My Novel NEPC Gene Set identified using Beltran, et al., over SB NEPC signature")
saveGSEAAnalysisOfNovelGeneSetAsPDF(sb_data)

printVIPERHeatmapMinusCIS(
  sb_data,
  g_mat,
  sig_list,
  cis_table.per_library,
  rna_seq_table,
  hallmarks_hm,
  kegg_hm,
  reactome_hm,
  filter_NP = TRUE
)
printVIPERHeatmapMinusCIS(
  sb_data,
  g_mat,
  sig_list,
  cis_table.per_library,
  rna_seq_table,
  hallmarks_hm,
  kegg_hm,
  reactome_hm,
  filter_NP = FALSE
)

