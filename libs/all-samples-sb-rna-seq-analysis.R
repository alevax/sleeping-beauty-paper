##
# Gene Expression Analysis on Pilot + New Cohort of SB RNA-Seq samples
# --------------------------------------------------------------------
# source("sources/sleeping-beauty/all-samples-sb-rna-seq-analysis.R")

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(AnnotationDbi))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(viper))
suppressPackageStartupMessages(library(crayon))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(org.Mm.eg.db))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(pheatmap))

#####################################################################################
#################################                   #################################
#################################    FILE PATHS     #################################
#################################                   ################################# 
#####################################################################################

filepath_SB_LIST_CIS_and_RNAseq_2020_OCT_20_xlsx = "data/mouse-analysis/SB List CIS and RNAseq_2020_OCT_20.xlsx"
filepath_sb_rna_seq_samples_table_updated_Mar2022_xlsx = "data/pathology/sb-rna-seq-samples-table_updated_Mar2022.xlsx"
filepath_raw_counts_all_mice_gene_symbols_kallisto_counts_rds = "data/sb-counts/raw-counts-all-mice-gene-symbols-kallisto-counts.rds"
filepath_cis_all_nr_prostate_0_001_csv = "data/cis-results/cis_all-nr-prostate-0.001.csv"
filepath_gemms_network_unPruned_rds = "data/sb-networks/gemms-network_unPruned.rds"

# source("../vaxtools/R/utils.R")
# source("../vaxtools/R/cross-species-utils.R")
source("libs/tools/utils.R")
source("libs/tools/cross_species_utils.R")
source("libs/tools/net-enrichment.R")
source("libs/tools/gemms-pathway-analysis-aREA.R")

# source("libs/tools/interactome_handler.R")

# create_workspace(run_dir = "sb-gene-expression-analysis-with-n1-lists-interactome")






# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
# -------------------------------------------------------------------------------------------
# -------------------------- LOAD DATASET FROM NEW COHORT FUNCTIONS -------------------------
# -------------------------------------------------------------------------------------------
# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&

# ------------------------------------- HELPER FUNCTIONS ------------------------------------
load_dataset_from_new_cohort_using_excel_file <- function(){
  x <- readxl::read_excel(filepath_SB_LIST_CIS_and_RNAseq_2020_OCT_20_xlsx , sheet = 5 , col_types = "text")
  colnames(x) <- make.names(colnames(x))
  x <- x[ grepl("CMZ" , x$RNAseq.ID ) , ]
  
  y <- readxl::read_excel(filepath_SB_LIST_CIS_and_RNAseq_2020_OCT_20_xlsx , sheet = 4 , col_types = "text")
  colnames(y) <- make.names(colnames(y))			
  y <- y[ grepl("CMZ" , y$RNAseq.ID ) , ]
  
  z <- inner_join( x , y )
  
  z$Bioanalyzer.ID <- NULL
  z$Quality.control <- NULL
  # z$Metastasis.sample <- NULL
  
  z$qPCR.cga <- z$qPCR.cga %>% as.numeric() %>% round(3)
  z$qPCR.Syn <- z$qPCR.Syn %>% as.numeric() %>% round(3)
  
  z <- z %>% dplyr::rename("phenotype"="NEPC.phenotype")
  z$phenotype <- ifelse(z$phenotype,"NEPC","Adeno")
  z <- z %>% dplyr::rename("castration_status"="Castration")
  z$castration_status <- ifelse(z$castration_status,"castrated","intact")	
  
  z$sb_status <- ifelse( grepl( "\\(\\+\\)" , z$Model ) , "activated" , "not-active" )
  z$genotype <- ifelse( grepl( "NPp53" , z$Model ) , "NPp53" , "NP" )			
  
  z <- z %>% dplyr::rename("has_dna_paired_samples"="CIS.analysis")
  z$has_dna_paired_samples <- ifelse(z$has_dna_paired_samples,"Yes","No")				
  sb_all_rna_metadata <- z %>% dplyr::rename("sample_id"="RNAseq.ID")			
  
  

  sb_all_rna_metadata_updatedpathology <-
    readxl::read_xlsx(filepath_sb_rna_seq_samples_table_updated_Mar2022_xlsx)
  stopifnot(sb_all_rna_metadata_updatedpathology$sample_id == sb_all_rna_metadata$sample_id)
  sb_all_rna_metadata_updatedpathology$phenotype <- sb_all_rna_metadata_updatedpathology$phenotype %>%
    dplyr::recode("NEPC (Group 2)" = "NEPC", "non-NEPC (Group 1)" = "Adeno")
  # View(data.frame("sample_id" = sb_all_rna_metadata$sample_id, "phenotype_old" = sb_all_rna_metadata$phenotype, "phenotype_new" = sb_all_rna_metadata_updatedpathology$phenotype))
  sb_all_rna_metadata$phenotype = sb_all_rna_metadata_updatedpathology$phenotype

  writexl::write_xlsx(sb_all_rna_metadata, "experiments/all-samples-sb-rna-seq-analysis/reports/sb-rna-seq-samples-table.xlsx")
  
  
  
  counts <- readRDS(filepath_raw_counts_all_mice_gene_symbols_kallisto_counts_rds)
  str(counts,1)
  
  index <- colnames(counts) %in% sb_all_rna_metadata$sample_id
  dim(sb_all_rna_metadata)
  counts <- counts[,index]
  dim(counts)
  
  index <- match( sb_all_rna_metadata$sample_id , colnames(counts) )
  counts <- counts[,index]
  stopifnot(identical(sb_all_rna_metadata$sample_id , colnames(counts)))	
  
  
  sb_data <- list()
  sb_data$raw_counts <- counts
  sb_data$metadata <- sb_all_rna_metadata
  return(sb_data)
}

# -------------------------------------- MAIN FUNCTIONS -------------------------------------
load_dataset_from_new_cohort <- function(){
  # sb_analysis.filename <- file.path("data/mouse-analysis/processed_data/sb-rna-seq-data-analysis.rds")
  # sb_analysis.filename <- file.path("data/oncomatch-analysis/SleepingBeauty/processed_data/sb-rna-seq-data-analysis.rds")
  # if(!file.exists(sb_analysis.filename)){
    sb_data <- load_dataset_from_new_cohort_using_excel_file()
    saveRDS(sb_data, "experiments/all-samples-sb-rna-seq-analysis/processed_data/sb-rna-seq-data-analysis.rds")
  # } else {
  #   sb_data <- readRDS(sb_analysis.filename)
  #   print_msg_info(">>> >> Loaded processed data file")
  #   print_msg_info(str(sb_data,1))
  # }
  return(sb_data)
}

# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
# -------------------------------------------------------------------------------------------
# ------------------------------------ PCA PLOT FUNCTIONS -----------------------------------
# -------------------------------------------------------------------------------------------
# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&

# ------------------------------------- HELPER FUNCTIONS ------------------------------------
getPCA1PCA2Variance <- function(vst, n_top = 500){
  x <- assay(vst)
  x <- x[ order( rowVars(x) , decreasing = TRUE)[1:n_top] , ]
  x <- prcomp( t(x) ) # , center = TRUE , scale. = TRUE )
  percentVar <- summary(x)$importance[2,]
  PC1_var <- round( percentVar["PC1"] * 100 / sum(percentVar) , 2 )
  PC2_var <- round( percentVar["PC2"] * 100 / sum(percentVar) , 2 )
  pca_var_list <- list(PC1_var, PC2_var)
  names(pca_var_list) <- c("PC1_var", "PC2_var")
  return(pca_var_list)
}
performPCAUsingDESeq2 <- function(sb_data, n_top = 500){
  require(DESeq2)
  require(ggrepel)
  dds <- DESeq2::DESeqDataSetFromMatrix( sb_data$raw_counts, colData = sb_data$metadata , design = ~1)
  vst <- DESeq2::vst( object = dds , blind = TRUE )
  
  pca.results <- DESeq2::plotPCA( vst ,
                                  ntop = n_top ,
                                  # ntop = nrow(assay(vst)) ,
                                  intgroup = c("phenotype","genotype","sb_status") , # ,"treatment_type","sampling_hour","purity") ,
                                  returnData = TRUE )
  pca_var_list <- getPCA1PCA2Variance(vst, n_top)
  
  pca_list <- c(list(pca.results), pca_var_list)
  names(pca_list)[1] <- "pca.results"
  return(pca_list)
}
getPCAPlot <- function(pca_list){
  pca.results <- pca_list[["pca.results"]]
  PC1_var <- pca_list[["PC1_var"]]
  PC2_var <-pca_list[["PC2_var"]]
  
  pca_plot <- ggplot( pca.results , aes( PC1 , PC2 , color = sb_status , shape = phenotype ) ) +
    # geom_text_repel( aes( label = sb_status ) , size = 3 ) +
    geom_point( alpha = 0.75 , size = 3 ) +
    xlab( paste0( "PC1: " , PC1_var , "% variance") ) +
    ylab( paste0( "PC2: " , PC2_var , "% variance") ) +
    coord_fixed() +
    theme_light() +
    facet_wrap(~genotype)
  # theme( legend.position = "none" )
  return(pca_plot)
}
getPCAPlot_NPp53 <- function(pca_list){
  pca.results <- pca_list[["pca.results"]] %>%
    filter(genotype == "NPp53")
  PC1_var <- pca_list[["PC1_var"]]
  PC2_var <-pca_list[["PC2_var"]]
  
  # dark orange
  # no grid
  # remove header
  # size of Rstudio panel - set width and height, 3x3, 4x4
  pca_plot_NPp53 <- ggplot( pca.results , aes( PC1 , PC2 , color =  phenotype, shape = sb_status ) ) +
    # geom_text_repel( aes( label = sb_status ) , size = 3 ) +
    geom_point( alpha = 0.75 , size = 3 ) +
    xlab( paste0( "PC1: " , PC1_var , "% variance") ) +
    ylab( paste0( "PC2: " , PC2_var , "% variance") ) +
    coord_fixed() +
    theme_light() +
    scale_shape_manual(values = c(17,19)) + 
    scale_color_manual(values = c("NEPC" = '#1f77b4', "Adeno" = '#ED820E')) + #'#ff7f0e'))
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA))
  # theme( legend.position = "none" )
  return(pca_plot_NPp53)
}
getPCAPlot_NP_NPp53 <- function(pca_list){
  pca.results <- pca_list[["pca.results"]]
  PC1_var <- pca_list[["PC1_var"]]
  PC2_var <-pca_list[["PC2_var"]]
  
  # shape by the genotype
  pca_plot_NP_NPp53 <- ggplot( pca.results , aes( PC1 , PC2 , color =  phenotype, shape = genotype ) ) +
    # geom_text_repel( aes( label = sb_status ) , size = 3 ) +
    geom_point( alpha = 0.75 , size = 3 ) +
    xlab( paste0( "PC1: " , PC1_var , "% variance") ) +
    ylab( paste0( "PC2: " , PC2_var , "% variance") ) +
    coord_fixed() +
    theme_light() +
    scale_shape_manual(values = c(19, 17)) + 
    scale_color_manual(values = c("NEPC" = '#1f77b4', "Adeno" = '#ED820E')) + #'#ff7f0e'))
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA))
  # theme( legend.position = "none" )
  return(pca_plot_NP_NPp53)
}

savePCAPlot <- function(pca_plot, n_top = 500){
  pdf( file = file.path( "experiments/all-samples-sb-rna-seq-analysis/reports/PCA/" ,
                         paste0("gene-expression-pca-analysis-ntopgenes-",n_top,".pdf") ) , width = 6 , height = 3 )
  print(pca_plot)
  dev.off()
  
  jpeg(file.path("experiments/all-samples-sb-rna-seq-analysis/reports/PCA/",
                 paste0("gene-expression-pca-analysis-ntopgenes-",n_top,".png")) , width = 6 , height = 3 , units = "in" , res = 300 )
  print(pca_plot)
  dev.off()
}

# -------------------------------------- MAIN FUNCTIONS -------------------------------------
plotPCAUsingDESeq2 <- function(sb_data){
  dir.create("experiments/all-samples-sb-rna-seq-analysis/reports/PCA/")
  stopifnot( identical( colnames(sb_data$raw_counts) , sb_data$metadata$sample_id ) )
  
  n_top <- 500
  pca_list <- performPCAUsingDESeq2(sb_data, n_top)
  pca_plot <- getPCAPlot(pca_list)
  
  pca_plot_NP_NPp53 <- getPCAPlot_NP_NPp53(pca_list)
  pdf( file = paste0( "experiments/all-samples-sb-rna-seq-analysis/reports/PCA/gene-expression-pca-analysis-ntopgenes-",n_top,"-NP_NPp53.pdf"), width = 4 , height = 4 )
  print(pca_plot_NP_NPp53)
  dev.off()
  
  pca_plot_NPp53 <- getPCAPlot_NPp53(pca_list)
  pdf( file = paste0( "experiments/all-samples-sb-rna-seq-analysis/reports/PCA/gene-expression-pca-analysis-ntopgenes-",n_top,"-NPp53.pdf"), width = 4 , height = 4 )
  print(pca_plot_NPp53)
  dev.off()
  
  message("- Printing PCA on PDF ...")
  savePCAPlot(pca_plot, n_top)
}	

# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
# -------------------------------------------------------------------------------------------
# ------------------------------- GATHER CIS ASSOCIATED GENES -------------------------------
# -------------------------------------------------------------------------------------------
# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
gatherCISAssociatedGenes <- function(cis_table){
  colnames(cis_table) <- make.names(colnames(cis_table))
  cis_table$cis_id <- paste0("CIS_",1:nrow(cis_table))
  
  ## Fisher's integration of 3 different pvalues estimation (this command has to come first for a proper Fisher's test)
  cis_table$fisher_pvalue <- apply( cbind( cis_table$pvalueinsert. , cis_table$pvaluelibrary. , cis_table$pvalueregion. ) , 1 , function(x) metap::sumlog(x)$p )
  ## Removing trailing , (commas)
  cis_table$gene_name <- cis_table$gene_name %>% map_chr( function(x) gsub( ",$", "" , x ) )
  ## Separating rows by making a row for each Gene and copying the rest of the columns
  cis_table <- cis_table %>% separate_rows(gene_name)
  ## Adding row for human orthologous gene (starting from mouse genes)
  cis_table$human_gene_name <- cis_table$gene_name %>% mouse_to_human(na.rm = FALSE)
  cis_table <- cis_table %>% dplyr::select(cis_id,fisher_pvalue,mouse_gene_name=gene_name,human_gene_name,library_name,everything())
  
  print_msg_info(">>> >> > ! Found " , 
                 length(unique( cis_table$mouse_gene_name[ !is.na(cis_table$mouse_gene_name) ] )) ,
                 " unique CIS-associated genes")
  return(cis_table)
}


# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
# -----------------------------------------------------------------------------------
# ----------------------- GENE EXPRESSION WITH EDGER FUNCTIONS ----------------------
# -----------------------------------------------------------------------------------
# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&

# --------------------------------- HELPER FUNCTION ---------------------------------
rankNorm <- function(x,FUN=median,trim=0)
{
  .rank <- apply(x, 2, rank)
  .median <- apply(.rank, 1, FUN, trim = trim )
  .mad <- apply(.rank, 1, mad)
  x <- (.rank - .median)/.mad
  
  message("- Number of NA features: " , sum( rowSums(is.na(x)) ) ) 
  message("- Number of Inf features: " , sum( rowSums(x) == Inf ) ) 
  message("- Number of 0 features: " , sum( rowSums(x) == 0 ) ) 
  
  message("- Features to Remove:")
  .null.features <- x[ which( rowSums( apply(x,2,is.na) ) > 0 ) , ]
  # print(.null.features)
  message("- Removing NULL/NA features ...")
  x <- na.omit(x)
  message("- Number of NA features: " , sum( rowSums(is.na(x)) ) )
  message("- Number of Inf features: " , sum( rowSums(x) == Inf ) )
  message("- Number of 0 features: " , sum( rowSums(x) == 0 ) ) 
  
  return(x)
}  
saveCountsMatricesToExcel <- function(cpm, raw_counts){
  dir.create("experiments/all-samples-sb-rna-seq-analysis/processed_data/counts/", showWarnings = FALSE)
  print_msg_info(">>> >> >> Saving CPM Matrix to Excel")
  cpm <- cpm( raw_counts[ rowSums(raw_counts) > 0 ,  ] , log = FALSE )	
  print(colSums(cpm)[1:5])
  cpm <- round(cpm,2)
  cpm.df <- as.data.frame(cpm)
  cpm.df$gene <- rownames(cpm)
  cpm.df <- cpm.df[ , rev(colnames(cpm.df)) ]
  writexl::write_xlsx( cpm.df , file.path("experiments/all-samples-sb-rna-seq-analysis/processed_data/counts/cpm-43-sb-mice.xlsx") )
  
  print_msg_info(">>> >> >> Saving scaled CPM Matrix to Excel")
  cpm <- cpm( raw_counts[ rowSums(raw_counts) > 0 ,  ] , log = FALSE )	
  print(colSums(cpm)[1:5])
  cpm <- t(scale(t(cpm)))
  cpm <- round(cpm,2)
  cpm.df <- as.data.frame(cpm)
  cpm.df$gene <- rownames(cpm)
  cpm.df <- cpm.df[ , rev(colnames(cpm.df)) ]
  writexl::write_xlsx( cpm.df , file.path("experiments/all-samples-sb-rna-seq-analysis/processed_data/counts/cpm-scaled-43-sb-mice.xlsx") )			
  
  print_msg_info(">>> >> >> Saving rankNorm CPM Matrix to Excel")
  cpm <- cpm( raw_counts[ rowSums(raw_counts) > 0 ,  ] , log = FALSE )	
  print(colSums(cpm)[1:5])
  cpm <- rankNorm(cpm)
  cpm <- round(cpm,2)
  cpm.df <- as.data.frame(cpm)
  cpm.df$gene <- rownames(cpm)
  cpm.df <- cpm.df[ , rev(colnames(cpm.df)) ]
  writexl::write_xlsx( cpm.df , file.path("experiments/all-samples-sb-rna-seq-analysis/processed_data/counts/cpm-ranknorm-43-sb-mice.xlsx") )
}
saveVolcanoPlotOfDGEGenes <- function(my_tibble, anEffect){
  require(ggplot2)
  require(tidyverse)
  require(ggrepel)
  
  my_plot <- ggplot( data = my_tibble , aes( x = logFC , y = logFDR ) ) +
    geom_point( color = "gray" ) +
    geom_point( alpha = 0.5 , stroke = 1 , size = 2 , color = "gray" , shape = 20 ) +
    geom_point( data = my_tibble %>% filter( is_diff_expr == "Yes" ) , color = "red" , alpha = 0.5 , stroke = 2 , size = 3 , fill = NA , shape = 1 )  +
    geom_label_repel( data = my_tibble %>% filter( is_diff_expr == "Yes" ) , mapping = aes( label = gene ) ) +
    
    xlab( paste0( "log2(Fold Change)") ) +
    ylab( paste0( "-log10(FDR)" ) ) +
    ggtitle( paste0( "Volcano Plot") ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size=15,face="bold") , # family="Open Sans") ,
      axis.title.x = element_text(size=15,face="bold") , # family="Open Sans") ,
      axis.title.y = element_text(size=15,face="bold") ) # family="Open Sans") )
  
  
  pdf(file.path(paste0("experiments/all-samples-sb-rna-seq-analysis/reports/edgeR/", anEffect, "/"),
                paste0(anEffect,"-volcano-plot-dge-genes.pdf")) )
  print(my_plot)
  dev.off()
  
  jpeg(file.path(paste0("experiments/all-samples-sb-rna-seq-analysis/reports/edgeR/", anEffect, "/"),
                 paste0(anEffect,"-volcano-plot-dge-genes.jpg")) , width = 6 , height = 6 , units = "in" , res = 300 )
  print(my_plot)
  dev.off()							
}
saveVolcanoPlotOfMarkerGenes <- function(my_tibble, anEffect){
  require(ggplot2)
  require(tidyverse)
  require(ggrepel)
  
  my_plot <- ggplot( data = my_tibble , aes( x = logFC , y = logFDR ) ) +
    geom_point( color = "gray" ) +
    geom_point( alpha = 0.5 , stroke = 1 , size = 2 , color = "gray" , shape = 20 ) +
    geom_point( data = my_tibble %>% filter( is_a_known_marker == "Yes" ) , color = "red" , alpha = 0.5 , stroke = 2 , size = 3 , fill = NA , shape = 1 )  +
    geom_label_repel( data = my_tibble %>% filter( is_a_known_marker == "Yes" ) , mapping = aes( label = gene ) ) +
    
    xlab( paste0( "log2(Fold Change)") ) +
    ylab( paste0( "-log10(FDR)" ) ) +
    ggtitle( paste0( "Volcano Plot") ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size=15,face="bold") , # family="Open Sans") ,
      axis.title.x = element_text(size=15,face="bold") , # family="Open Sans") ,
      axis.title.y = element_text(size=15,face="bold") ) # family="Open Sans") )
  
  pdf(file.path(paste0("experiments/all-samples-sb-rna-seq-analysis/reports/edgeR/", anEffect, "/"),
                paste0(anEffect,"-volcano-plot-marker-genes.pdf")) )
  print(my_plot)
  dev.off()
  
  jpeg(file.path(paste0("experiments/all-samples-sb-rna-seq-analysis/reports/edgeR/", anEffect, "/"),
                 paste0(anEffect,"-volcano-plot-marker-genes.jpg")) , width = 6 , height = 6 , units = "in" , res = 300 )
  print(my_plot)
  dev.off()								
}
saveVolcanoPlotOfCISGenes <- function(my_tibble, anEffect){
  require(ggplot2)
  require(tidyverse)
  require(ggrepel)
  
  my_plot <- ggplot( data = my_tibble , aes( x = logFC , y = logFDR ) ) +
    geom_point( color = "gray" ) +
    geom_point( alpha = 0.5 , stroke = 1 , size = 2 , color = "gray" , shape = 20 ) +
    geom_point( data = my_tibble %>% filter( is_cis_gene == "Yes" & FDR < 1e-02 & abs(logFC) > 1 ) , color = "red" , alpha = 0.5 , stroke = 2 , size = 3 , fill = NA , shape = 1 )  +
    geom_label_repel( data = my_tibble %>% filter( is_cis_gene == "Yes" & FDR < 1e-02 & abs(logFC) > 1 ) , mapping = aes( label = gene ) , force = 10 ) +
    
    xlab( paste0( "log2(Fold Change)") ) +
    ylab( paste0( "-log10(FDR)" ) ) +
    ggtitle( paste0( "Volcano Plot") ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size=15,face="bold") , # family="Open Sans") ,
      axis.title.x = element_text(size=15,face="bold") , # family="Open Sans") ,
      axis.title.y = element_text(size=15,face="bold") ) # family="Open Sans") )
  
  pdf(file.path(paste0("experiments/all-samples-sb-rna-seq-analysis/reports/edgeR/", anEffect, "/"),
                paste0(anEffect,"-volcano-plot-cis-genes.pdf")) )
  print(my_plot)
  dev.off()
  
  jpeg(file.path(paste0("experiments/all-samples-sb-rna-seq-analysis/reports/edgeR/", anEffect, "/"),
                 paste0(anEffect,"-volcano-plot-cis-genes.jpg")) , width = 6 , height = 6 , units = "in" , res = 300 )
  print(my_plot)
  dev.off()				
}
lapplyToMyContrasts <- function(contrast.matrix, fit){
  dir.create("experiments/all-samples-sb-rna-seq-analysis/reports/edgeR/")
  my_contrasts <- colnames(contrast.matrix)
  res <- lapply( my_contrasts , function(anEffect) {
    
    qlf <- glmQLFTest(fit, contrast = contrast.matrix[,anEffect] )
    tt <- topTags( qlf , sort.by = "PValue" , n = nrow(qlf) )
    ges <- qnorm( tt$table$PValue/2 , lower.tail = FALSE ) * sign(tt$table$logFC)
    names(ges) <- rownames(tt)
    
    de <- decideTestsDGE( qlf , adjust.method = "BH", p.value = 0.1 )
    up_genes <-  rownames(de)[ de == +1 ]
    down_genes <-  rownames(de)[ de == -1 ]
    # plotSmear(x$lrt, de.tags=c(up_genes,down_genes))
    # abline(h = c(-2, 2), col = "blue")
    print_msg_info(anEffect)
    print_msg_info(table(de))
    
    my_tibble <- qlf$table %>% rownames_to_column("gene") %>% as_tibble()
    my_tibble$gene <- rownames(qlf$table)
    my_tibble$FDR <- p.adjust( my_tibble$PValue , method = "BH" )
    my_tibble$logFDR <- -log10( my_tibble$FDR )
    
    # my_tibble %>% arrange(-FDR) %>% slice_tail(n = 20)
    x1 <- my_tibble %>% arrange(logFC) %>% slice_tail(n = 10) %>% dplyr::select(gene) %>% pull()
    x2 <- my_tibble %>% arrange(logFDR) %>% slice_tail(n = 10) %>% dplyr::select(gene) %>% pull()
    top_dge_genes <- unique(c(x1,x2))
    
    my_tibble$is_diff_expr <- ifelse( my_tibble$gene %in% top_dge_genes , "Yes" , "No" )
    my_tibble$is_cis_gene <- ifelse( my_tibble$gene %in% cis_genes , "Yes" , "No" )
    my_tibble$is_a_known_marker <- ifelse( my_tibble$gene %in% known_nepc_markers , "Yes" , "No" )
    
    dir.create(paste0("experiments/all-samples-sb-rna-seq-analysis/reports/edgeR/", anEffect, "/"))
    saveRDS( my_tibble , file.path("experiments/all-samples-sb-rna-seq-analysis/reports/edgeR/",
                                   paste0(anEffect, "/", anEffect,"-tibble.rds")) )
    
    print_msg_info(">>> Volcano Plot of DGE Genes")
    saveVolcanoPlotOfDGEGenes(my_tibble, anEffect)
    
    print_msg_info(">>> Volcano Plot of Marker Genes")
    saveVolcanoPlotOfMarkerGenes(my_tibble, anEffect)
    
    print_msg_info(">>> Volcano Plot of CIS Genes")
    saveVolcanoPlotOfCISGenes(my_tibble, anEffect)
    .ret <- list( anEffect , qlf , tt , ges )
    names(.ret) <- c("anEffect","qlf","tt","ges")
    return(.ret)
    
  })
  names(res) <- my_contrasts
  return(res)
}

# ---------------------------------- MAIN FUNCTION ----------------------------------
geneExpressionAnalysisWithEdgeR <- function(sb_data){
  raw_counts <- sb_data$raw_counts
  sb_metadata <- sb_data$metadata
  
  require(edgeR)
  cpm_log <- cpm( raw_counts , log = TRUE )
  print_msg_info(">>> >> Saving counts Matrices to Excel")
  saveCountsMatricesToExcel(cpm, raw_counts)
  
  median_log2_cpm <- apply(cpm_log, 1, median)
  hist(median_log2_cpm)
  expr_cutoff <- -3
  abline(v = expr_cutoff, col = "red", lwd = 3)
  sum(median_log2_cpm > expr_cutoff)
  raw_counts_filtered <- raw_counts[median_log2_cpm > expr_cutoff, ]
  dim(raw_counts_filtered)
  
  set.seed(666)
  
  sb_metadata$groups <- paste( sb_metadata$genotype , sb_metadata$sb_status , sb_metadata$phenotype , sep = "_" )
  sb_metadata$groups <- make.names(sb_metadata$groups)
  sb_metadata$groups <- factor(sb_metadata$groups)
  sb_metadata$groups <- relevel(sb_metadata$groups,ref="NP_not.active_Adeno")
  
  y <- DGEList( counts = raw_counts_filtered , group = sb_metadata$groups )
  y <- calcNormFactors( y , method = "RLE")
  # y <- calcNormFactors( y , method = "TMM" , refColumn = sb_metadata$groups == "NP_not.active_Adeno")
  
  colors <- brewer.pal( nlevels(sb_metadata$groups) , "Set1")
  names(colors) <- as.character(levels(y$samples$group))
  plotMDS(y, method = "bcv", col = colors[ y$samples$group ] )
  legend("bottom", names(colors) , col = colors , pch=20 )
  
  design <- model.matrix( ~0 + groups , data = sb_metadata ) # + batch + rin)
  colnames(design) <- gsub("groups","",colnames(design))
  colnames(design)
  
  y <- estimateGLMTrendedDisp(y, design = design)
  y <- estimateGLMTagwiseDisp(y, design = design)
  
  
  
  contrast.matrix <- makeContrasts( levels = colnames(design) ,
                                    NP_SB_effect_on_Adeno = NP_activated_Adeno - NP_not.active_Adeno ,
                                    NPp53_SB_effect_on_Adeno = NPp53_activated_Adeno - NPp53_not.active_Adeno ,
                                    NPp53_SB_effect_on_NEPC = NPp53_activated_NEPC - NPp53_not.active_Adeno ,
                                    NPp53_NEPC_effect_on_SB_activated = NPp53_activated_NEPC - NPp53_activated_Adeno
  )
  
  fit <- glmQLFit(y, design)
  res <- lapplyToMyContrasts(contrast.matrix, fit)
  
  # str(res,2)
  
  # genes <- sort(rownames(cpm(y)))
  genes <- sort( unique( unlist( lapply( res , function(x) names(x$ges) ) ) ) )
  gesmat <- lapply( lapply( res , `[[` , "ges" ) , function(x) x[ genes ] )
  # str(gesmat)
  gesmat <- do.call( cbind , gesmat )		
  return(gesmat)
}


# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
# -----------------------------------------------------------------------------------
# ---------------------------- GEXPR COMPLEXHEATMAP FUNCS ---------------------------
# -----------------------------------------------------------------------------------
# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
saveGeneExpressionHeatmap <- function(sb_data){
  require(circlize)
  require(ComplexHeatmap)
  df_annot.tibble <- sb_data$metadata
  mat <- cpm( sb_data$raw_counts , log = FALSE )
  stopifnot( identical( colnames(mat) , df_annot.tibble$sample_id ) )
  
  df_annot.tibble <- df_annot.tibble %>% dplyr::select(genotype,sb_status,phenotype)
  
  # df_annot.tibble <- data.frame( HISTOLOGY = df_annot.tibble$nepc_phenotype ,
  # 															 SB_STATUS = df_annot.tibble$sb_status
  # 														# SYP = df_annot.tibble$qPCR.Syn %>% as.numeric() ,
  #                            # CGA = df_annot.tibble$qPCR.Cga %>% as.numeric() ,
  #                             )
  
  col_fun = colorRamp2(c(-3, 0, 3), c("deepskyblue3", "white", "brown3"))
  
  df_annot.df <- as.data.frame(df_annot.tibble)
  
  df_annot_cols = HeatmapAnnotation( df = df_annot.df ,
                                     col = list( phenotype = c("NEPC"="steelblue","Adeno"="tan1") ,
                                                 genotype = c("NP"="green4","NPp53"="purple") ,
                                                 sb_status = c("activated"="red3","not-active"="gray50")
                                                 # SYP = colorRamp2(c(0,max(df_annot.df$SYP,na.rm=T) ), c("white", "darkblue")) ,
                                                 # CGA = colorRamp2(c(0,max(df_annot.df$CGA,na.rm=T) ), c("white", "darkblue")) 
                                     )
  )  
  
  index <- order(rowVars(mat),decreasing = T)[1:50]
  mat <- mat[index,]
  
  mat <- t(scale(t(mat)))
  
  main_hm <- Heatmap( mat ,
                      # col = circlize::colorRamp2(c(-3, 0, 3), c( muted("blue"), "white", muted("red"))),
                      col = col_fun ,
                      heatmap_legend_param = list(color_bar = "continuous") ,
                      name = "Gene Expression" )
  
  ht_list = df_annot_cols %v% main_hm
  # ht_list = df_annot_cols %v% df_annot_cols_bottom %v% main_hm
  p <- draw(ht_list, column_km = 1 )
  
  pdf( file.path( "experiments/all-samples-sb-rna-seq-analysis/reports/gene-expression-heatmap.pdf") ) # , width = 6, height = 6 )
  print(p)
  dev.off()  		
}

# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
# -----------------------------------------------------------------------------------
# ------------------------------- VIPER ANALYSIS FUNCS ------------------------------
# -----------------------------------------------------------------------------------
# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
saveLargeMatToExcel <- function(mat, filename){
  mat_xlsx <- data.frame(cbind(rownames(mat), mat))
  colnames(mat_xlsx)[1] <- "rowname"
  writexl::write_xlsx(mat_xlsx, filename, col_names = TRUE)
}
addVIPERAnalysisToSBData <- function(raw_counts, mouse_interactome, use_NP_as_reference = FALSE){
  require(edgeR)
  cpm_log <- cpm( raw_counts , log = TRUE )
  
  # sb_vpmat <- viper(cpm_log,pruneRegulon(mouse_interactome,50),method = "mad")
  # sb_data$vpmat <- sb_vpmat
  # attr(sb_data$vpmat,"pruneRegulon") <- "50"
  # attr(sb_data$vpmat,"signature_method") <- "viper_mad"				
  
  ss_gesmat <- rankNorm(cpm_log)
  saveLargeMatToExcel(ss_gesmat, "experiments/all-samples-sb-rna-seq-analysis/reports/sb_gesmat.xlsx")
  
  if(use_NP_as_reference){
    ss_vpsig <- viperSignature(eset = ss_gesmat, ref = ss_gesmat[, sb_data$metadata$genotype=="NP"])
    sb_vpmat <- viper( ss_vpsig$signature , pruneRegulon(mouse_interactome,50) , method = "none" )
  } else {
    sb_vpmat <- viper( ss_gesmat , pruneRegulon(mouse_interactome,50) , method = "none" )
  }
  saveLargeMatToExcel(sb_vpmat, "experiments/all-samples-sb-rna-seq-analysis/reports/sb_vpmat.xlsx")
  
  sb_data$vpmat <- sb_vpmat
  attr(sb_data$vpmat,"pruneRegulon") <- "50"
  attr(sb_data$vpmat,"signature_method") <- "rankNorm"
  return(sb_data)
}

# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
# -----------------------------------------------------------------------------------
# ------------------------------ PATHWAY ANALYSIS FUNCS -----------------------------
# -----------------------------------------------------------------------------------
# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
addPathwayAnalysisToSBData <- function(sb_data){
  # source("../vaxtools/R/pathway-analysis.R")
  # sb_pathways <- do_pathways_with_aREA(sb_vpmat %>% mouse_to_human() , threshold = 1 )
  pathways_dir <- "experiments/all-samples-sb-rna-seq-analysis/reports/pathways/"
  dir.create(pathways_dir)
  sb_pathways <- do_pathways_with_aREA(sb_data$vpmat %>% mouse_to_human(), out_dir = pathways_dir)
  sb_data$pathways <- sb_pathways
  sb_data$metadata$beltran_enrichment <- getNETenrichment( sb_data$vpmat %>% mouse_to_human() , 
                                                           regulonObject = beltran.net.regulon.obj )
  return(sb_data)
}	

# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
# -----------------------------------------------------------------------------------
# ------------------------------ GET HEATMAP LIST FUNCS -----------------------------
# -----------------------------------------------------------------------------------
# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&

# --------------------------------- HELPER FUNCTION ---------------------------------

getVIPERMatrixFromMarina <- function(gesmat, gemms_interactome){
  known_nepc_markers <- c("Ezh2","Ascl1","Onecut2","Foxa2","Foxa1","Sox2","Ar","Mycn","Sox11")
  
  
  for ( anEffect in colnames(gesmat) )
  {
    dir.create(paste0("experiments/all-samples-sb-rna-seq-analysis/reports/marina/", anEffect, "/"), showWarnings = FALSE)
    ms <- msviper(gesmat[,anEffect], pruneRegulon(gemms_interactome, 100) ) # PRUNE 100
    x <- sort(ms$es$nes)
    selected_mrs <- c(names(head(x,15)),names(tail(x,15))) %>% rev()
    pdf( file.path(paste0("experiments/all-samples-sb-rna-seq-analysis/reports/marina/", anEffect, "/"),
                   paste0(anEffect,"-signature-marina-analysis-mouse.pdf")) )
    plot(ms,selected_mrs)
    dev.off()
    
    pdf( file.path(paste0("experiments/all-samples-sb-rna-seq-analysis/reports/marina/", anEffect, "/"),
                   paste0(anEffect,"NEPC-signature-marina-analysis-mouse-known-nepc-markers.pdf")) , width = 4 )
    plot( ms , names(x[ names(x) %in% known_nepc_markers ]) %>% rev() )
    dev.off()			
  }
  
  # vpmat <- viper( gesmat , pruneRegulon(gemms_interactome,50) , method = "none")
  vpmat <- viper( gesmat , pruneRegulon(gemms_interactome, 100) , method = "none") # PRUNE 100
  return(vpmat)
}

# ---------------------------------- MAIN FUNCTION ----------------------------------

addMarinaAnalysisToSBData <- function(sb_data, gesmat, gemms_interactome){
  require(viper)
  dir.create("experiments/all-samples-sb-rna-seq-analysis/reports/marina/", showWarnings = FALSE)
  vpmat <- getVIPERMatrixFromMarina(gesmat, gemms_interactome)
  saveRDS(vpmat, file.path( "experiments/all-samples-sb-rna-seq-analysis/reports/vpmat.rds") )
  writexl::write_xlsx(data.frame(vpmat), "experiments/all-samples-sb-rna-seq-analysis/reports/vpmat_from_edgeR_signatures.xlsx")
  sb_data$vpmat_from_edgeR_signatures <- vpmat
  return(sb_data)
}

# # &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
# # -----------------------------------------------------------------------------------
# # ---------------------------- SAVE VPMAT USING NP AS REF ---------------------------
# # -----------------------------------------------------------------------------------
# # &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
# 
# saveVPMatWithNPSamplesAsReference <- function(sb_data){
#   cpm_log <- cpm( sb_data$raw_counts , log = TRUE )
#   ss_gesmat <- rankNorm(cpm_log)
#   ss_vpsig <- viperSignature(eset = ss_gesmat, ref = ss_gesmat[, sb_data$metadata$genotype=="NP"])
#   sb_vpmat_vpsig <- viper( ss_vpsig$signature , pruneRegulon(mouse_interactome,50) , method = "none" )
#   saveRDS(sb_vpmat_vpsig, "experiments/all-samples-sb-rna-seq-analysis/sb_vpmat_vpsig_npRef.rds")
# }


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

dir.create("experiments/all-samples-sb-rna-seq-analysis/")
dir.create("experiments/all-samples-sb-rna-seq-analysis/reports/")
dir.create("experiments/all-samples-sb-rna-seq-analysis/processed_data/")

## Loading Dataset From New Cohort ----
print_msg_info(">>> Loading Dataset From New Cohort")
sb_data <- load_dataset_from_new_cohort()

## PCA Plot Using DESeq2 ----
message(">>> PCA Plot Using DESeq2")
plotPCAUsingDESeq2(sb_data)

known_nepc_markers <- c("Ezh2","Ascl1","Onecut2","Foxa2","Foxa1","Sox2","Ar","Mycn","Sox11","Foxm1","Cenpf","Syp","Chga")

## Gathering CIS-assiciated genes ----
print_msg_info(">>> Gathering CIS-assiciated genes")
print_msg_warn(">>> *** This List is not derived from NE-only samples ***")
cis_table <- read_csv(filepath_cis_all_nr_prostate_0_001_csv, skip = 1)
cis_table <- gatherCISAssociatedGenes(cis_table)
cis_genes <- cis_table$mouse_gene_name
sb_data$cis_table <- cis_table
saveRDS(sb_data, "experiments/all-samples-sb-rna-seq-analysis/processed_data/sb-rna-seq-data-analysis.rds")

## Gene Expression Analysis with EdgeR ----
print_msg_info(">>> Gene Expression Analysis with EdgeR")
gesmat <- geneExpressionAnalysisWithEdgeR(sb_data)

## Heatmap at Gene Expression ----
print_msg_info(">>> Heatmap at Gene Expression Analysis with ComplexHeatmap")
saveGeneExpressionHeatmap(sb_data)

## VIPER Analysis ----
gemms_interactome <- readRDS(filepath_gemms_network_unPruned_rds)
class(gemms_interactome) <- "regulon"
print_msg_info(">>> Master Regulator Analysis with VIPER")
sb_data <- addVIPERAnalysisToSBData(sb_data$raw_counts, gemms_interactome, use_NP_as_reference = FALSE)

## Pathway Analysis ----
print_msg_info(">>> Pathway Analysis for GEMMs")
sb_data <- addPathwayAnalysisToSBData(sb_data)

## Master Regulator Analysis with Marina ----
print_msg_info(">>> Master Regulator Analysis with Marina")
sb_data <- addMarinaAnalysisToSBData(sb_data, gesmat, gemms_interactome)

## Saving SB Object ----
print_msg_info(">>> Saving SB Ojbect to: " , "experiments/all-samples-sb-rna-seq-analysis/processed_data/sb-rna-seq-data-analysis.rds" )
saveRDS(sb_data, "experiments/all-samples-sb-rna-seq-analysis/processed_data/sb-rna-seq-data-analysis.rds")

# saveVPMatWithNPSamplesAsReference(sb_data)
# 
