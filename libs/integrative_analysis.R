##
# SB Integrative Analysis for CIS / VIPER / CINDy
# -----------------------------------------------
# system.time({source("sources/sleeping-beauty/integrative-analysis.R")})

# Libraries
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(pryr))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(circlize))

# source("../vaxtools/R/utils.R")
suppressPackageStartupMessages(require(limma))
suppressPackageStartupMessages(require(matrixStats))
suppressPackageStartupMessages(require(crayon))
suppressPackageStartupMessages(require(tidyverse))

suppressPackageStartupMessages(library(AnnotationDbi))
suppressPackageStartupMessages(library(stats4))
suppressPackageStartupMessages(library(BiocGenerics))
suppressPackageStartupMessages(library(multtest))

suppressPackageStartupMessages(library(metap))


# source("../vaxtools/R/cross-species-utils.R")
# create_workspace("cis-integrative-analysis-usingNES_9_asThresholdForCINDyMRs-CINDyMergedSU2C+TCGA-gexFDR1em3")

#####################################################################################
#################################                   #################################
#################################    FILE PATHS     #################################
#################################                   ################################# 
#####################################################################################
filepath_cindy_table_su2c_plus_tcga_merged_rds <- "data/cindy-results/cindy-table-su2c-plus-tcga-merged.rds"
filepath_NPp53_NEPC_effect_on_SB_activated_tibble_rds <- "experiments/all-samples-sb-rna-seq-analysis/reports/edgeR/NPp53_NEPC_effect_on_SB_activated/NPp53_NEPC_effect_on_SB_activated-tibble.rds"
filepath_vpmat_rds <- "experiments/all-samples-sb-rna-seq-analysis/reports/vpmat.rds"
filepath_cis_all_nr_prostate_0_001_csv <- "data/cis-results/cis_all-nr-prostate-0.001.csv"
filepath_SB_List_CIS_and_RNAseq_2020_OCT_20_xlsx <- "data/mouse-analysis/SB List CIS and RNAseq_2020_OCT_20.xlsx"
filepath_tfs_csv <- "data/regulators-lists/tfs.csv"
filepath_cotfs_csv <- "data/regulators-lists/cotfs.csv"
filepath_sig_csv <- "data/regulators-lists/sig.csv"
filepath_ccle_tpm_symbols_rds <- "data/ccle-data/ccle-tpm-symbols.rds"
filepath_Cell_lines_annotations_20181226_txt <- "data/ccle-data/Cell_lines_annotations_20181226.txt"

source("libs/tools/cross_species_utils.R")
source("libs/tools/utils.R")

# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
# -------------------------------------------------------------------------------------------
# --------------------------------- PREPROCESSING FUNCTIONS ---------------------------------
# -------------------------------------------------------------------------------------------
# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
preprocessCindyTable <- function(cindy_table){
  cindy_table <- cindy_table %>% ungroup()
  cindy_table$Modulator.mouse <- cindy_table$Modulator %>% human_to_mouse(na.rm = F)
  
  x <- unique(cindy_table$significantTriplets)
  y <- rank(x)/max(rank(x))
  names(x) <- round((1-as.numeric(y))*99)+1
  cindy_table$ranking <- names( x[ match( cindy_table$significantTriplets , x ) ] )
  cindy_table$ranking <- paste0(cindy_table$ranking,"%")
  return(cindy_table)
}
preprocessMyTibble <- function(my_tibble, vpmat){
  # Overwriting the selected DGE for Volcano plot (we want all the DGE)
  my_tibble$is_diff_expr <- ifelse( my_tibble$FDR < 1e-3 , "Yes" , "No" )	
  
  x <- vpmat[,"NPp53_NEPC_effect_on_SB_activated"]
  my_tibble$NEPC_viper_score <- x[ match( my_tibble$gene , names(x) ) ]		
  
  mrs_candidates <- names(vpmat[,"NPp53_NEPC_effect_on_SB_activated"])[ vpmat[,"NPp53_NEPC_effect_on_SB_activated"] > 2 ]
  mrs_candidates.human <- mrs_candidates %>% mouse_to_human()	# Never used below	
  return(my_tibble)
}
preprocessCisTable <- function(cis_table){
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
  
  ## Hacking Gene Name because of bad mouse ---> human conversion ----
  cis_table$human_gene_name <- ifelse(cis_table$gene_name == "Mll3","KMT2C",cis_table$human_gene_name)
  cis_table$gene_name <- ifelse(cis_table$gene_name == "Mll3","Kmt2c",cis_table$gene_name)
  
  cis_table <- cis_table %>% dplyr::select(cis_id,fisher_pvalue,mouse_gene_name=gene_name,human_gene_name,library_name,everything())
  
  print_msg_info(">>> >> > ! Found " , 
                 length(unique( cis_table$mouse_gene_name[ !is.na(cis_table$mouse_gene_name) ] )) ,
                 " unique CIS-associated genes")
  return(cis_table)
}
preprocessSbExcel <- function(sb_excel){
  colnames(sb_excel) <- make.names(colnames(sb_excel))
  print_msg_info(">>> >> Found a total of " , 
                 length(unique(sb_excel$DNA.Sample.ID)) , 
                 " different DNA Barcodes")	
  
  print_msg_info(">>> >> Found a total of " , 
                 length(unique(cis_table$library_name[ cis_table$library_name %in% unique(sb_excel$DNA.Sample.ID) ])) ,
                 " different DNA Barcodes matched by RNA-Seq samples")	
  
  sb_excel$mouse_cohort <- ifelse( sb_excel$Pilot.cohort , "pilot" , "latest" )
  sb_excel$castration_status <- ifelse( sb_excel$Castration , "castrated" , "intact" )
  sb_excel$tumor_location <- ifelse( sb_excel$Metastasis.sample , "mets" , "primary" )
  sb_excel$phenotype <- ifelse( sb_excel$NE.Phenotype , "nepc" , "adeno" )
  return(sb_excel)
}
preprocessCCLE <- function(ccle.tpm, ccle_annotations){
  object_size(ccle.tpm)
  
  
  prostate_cell_lines.names <- ccle_annotations$CCLE_ID[ grepl( "PRAD" , ccle_annotations$tcga_code ) ]
  prostate_cell_lines.index <- colnames(ccle.tpm) %in% prostate_cell_lines.names
  colnames(ccle.tpm)[prostate_cell_lines.index]
  
  ccle_prad.tpm <- ccle.tpm[,prostate_cell_lines.index]
  dim(ccle_prad.tpm)
  
  ccle.log2tpm <- log2(ccle.tpm+1)
  ccle_prad.log2tpm <- log2(ccle_prad.tpm+1)
  return(ccle_prad.log2tpm)
}
# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
# -------------------------------------------------------------------------------------------
# ------------------------------- CIS TABLE PER LIBRARY FUNCS -------------------------------
# -------------------------------------------------------------------------------------------
# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&

# ------------------------------------- HELPER FUNCTIONS ------------------------------------
reformatCisTableWithEachRowPerLibrary <- function(cis_table.per_library){
  ## Separating rows by making a row for each Library and copying the rest of the columns
  print_msg_info(">>> >> > Separating rows by making a row for each Library and copying the rest of the columns")
  cis_table.per_library$library_name <- cis_table.per_library$library_name %>% map_chr( function(x) gsub( "::", "," , x ) )
  cis_table.per_library$library_name <- cis_table.per_library$library_name %>% map_chr( function(x) gsub( "^,|,$", "" , x ) )
  cis_table.per_library <- cis_table.per_library %>% separate_rows(library_name,sep = ",")
  
  cis_table.per_library <- cis_table.per_library %>%
    dplyr::select(cis_id,fisher_pvalue,mouse_gene_name,human_gene_name,library_name,everything())
  cis_table.per_library <- cis_table.per_library %>% arrange(fisher_pvalue)
  return(cis_table.per_library)
}
incorporateSbExcelDataIntoCisTablePerLibrary <- function(cis_table.per_library, sb_excel){
  cis_table.per_library$mouse_cohort <- sb_excel$mouse_cohort[ match( cis_table.per_library$library_name , sb_excel$DNA.Sample.ID ) ]
  cis_table.per_library$castration_status <- sb_excel$castration_status[ match( cis_table.per_library$library_name , sb_excel$DNA.Sample.ID ) ]
  cis_table.per_library$tumor_location <- sb_excel$tumor_location[ match( cis_table.per_library$library_name , sb_excel$DNA.Sample.ID ) ]
  cis_table.per_library$phenotype <- sb_excel$phenotype[ match( cis_table.per_library$library_name , sb_excel$DNA.Sample.ID ) ]
  
  cis_table.per_library <- cis_table.per_library %>% dplyr::select(cis_id,fisher_pvalue,mouse_gene_name,human_gene_name,phenotype,tumor_location,castration_status,mouse_cohort,library_name,everything())
  return(cis_table.per_library)
}

# -------------------------------------- MAIN FUNCTIONS -------------------------------------
saveCisTablePerLibrary <- function(cis_table, sb_excel){
  cis_table.per_library <- reformatCisTableWithEachRowPerLibrary(cis_table)
  cis_table.per_library <- incorporateSbExcelDataIntoCisTablePerLibrary(cis_table.per_library, sb_excel)
  my_filename <- file.path("experiments/integrative-analysis/processed_data/cis-table-per-library.rds")
  saveRDS(cis_table.per_library, my_filename)
}


# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
# -------------------------------------------------------------------------------------------
# ---------------------------------- DATA INTEGRATION FUNCS ---------------------------------
# -------------------------------------------------------------------------------------------
# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
integrateCindyTableAndCisTableIntoMyTibble <- function(my_tibble, cindy_table, cis_table){
  print_msg_info(">>> >> Integrating DATA")
  my_tibble$cis_fisher_pvalue <- cis_table$fisher_pvalue[ match( my_tibble$gene , cis_table$mouse_gene_name ) ]
  
  print_msg_warn(">>> ** CIS-associated Genes Excluded from RNA-Seq Analysis")
  # setdiff( cis_table$mouse_gene_name , my_tibble$gene_name )
  # length(setdiff( cis_table$mouse_gene_name , my_tibble$gene_name ))
  
  nes_threshold <- 9
  up <- my_tibble %>%
    filter(gene %in% human_to_mouse(cindy_table$TF)) %>%
    slice_max(NEPC_viper_score,n=25) %>%
    dplyr::select(gene) %>%
    pull()
  down <-  my_tibble %>%
    filter(gene %in% human_to_mouse(cindy_table$TF)) %>%
    slice_min(NEPC_viper_score,n=25) %>%
    dplyr::select(gene) %>%
    pull()
  top_and_bottom_cmrs <- c(up,down)
  writexl::write_xlsx(my_tibble %>%
                        filter(gene %in% top_and_bottom_cmrs) %>%
                        dplyr::select(proteins=gene,NES=NEPC_viper_score),
                      file.path("experiments/integrative-analysis/reports/nepc-top-50-candidate-mrs-table.xlsx") )	
  
  cindy_table.filtered <- cindy_table %>%
    filter(pvalue < 0.05) %>%
    filter(Modulator %in% cis_table$human_gene_name) %>%
    filter(TF %in% mouse_to_human(my_tibble$gene[abs(my_tibble$NEPC_viper_score) > nes_threshold]))
  
  cindy_table.filtered <- cindy_table.filtered %>%
    group_by(Modulator) %>%
    dplyr::mutate(significantTriplets_mean=mean(significantTriplets))
  cindy_table.filtered <- cindy_table.filtered %>%
    group_by(Modulator) %>%
    dplyr::mutate(cindy_pvalue_integrated=metap::sumlog(pvalue)[["p"]])
  # neg log pvalue addition due to overflow of p-values to 0
  cindy_table.filtered <- cindy_table.filtered %>%
    group_by(Modulator) %>%
    dplyr::mutate(cindy_neg_log_pvalue_integrated=metap::sumlog(pvalue, log.p = TRUE)[["p"]]*-1)
  
  cindy_table.filtered$TF.mouse <- human_to_mouse(cindy_table.filtered$TF, na.rm = FALSE)
  cindy_table.filtered <- left_join( cindy_table.filtered , 
                                     my_tibble %>% dplyr::select(gene,NEPC_viper_score) , by = c("TF.mouse"="gene") )
  writexl::write_xlsx(cindy_table.filtered,
                      file.path("experiments/integrative-analysis/reports/cindy-table-modulators-vs-tfs.xlsx") )
  
  print_msg_warn(">>> ** Using NES > abs(" , nes_threshold , ") for TFs as downstream effectors of Modulators")		
  print_msg_warn(">>> ** Using " , length(unique(cindy_table.filtered$TF)) , " TFs as downstream effectors of Modulators")		
  
  my_tibble <- left_join(my_tibble , 
                         cindy_table.filtered %>% dplyr::distinct(
                           Modulator.mouse,
                           cindy_pvalue_integrated,
                           cindy_neg_log_pvalue_integrated
                           ), by = c("gene"="Modulator.mouse") )
  
  x <- unique(cindy_table.filtered$Modulator) %>% human_to_mouse()
  my_tibble$is_candidate_modulator <- ifelse( is.na( x[ match( my_tibble$gene , x ) ] ) , "No" , "Yes" )
  my_tibble$modulator_ranking <- cindy_table.filtered$ranking[ match( my_tibble$gene , cindy_table.filtered$Modulator.mouse ) ]
  
  my_tibble$gene_human <- my_tibble$gene %>% mouse_to_human(na.rm = FALSE)
  my_tibble$is_c_mrs <- my_tibble$gene_human %in% unique(cindy_table.filtered$TF)
  sum(my_tibble$is_c_mrs)
  return(my_tibble)
}
# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
# -------------------------------------------------------------------------------------------
# ---------------------------------- NEPC SIGNATURES FUNCS ----------------------------------
# -------------------------------------------------------------------------------------------
# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
reformatMyTibbleToSaveNEPCSignatures <- function(my_tibble){
  my_tibble <- my_tibble %>% dplyr::select(gene,cis_fisher_pvalue,is_cis_gene,is_diff_expr,
                                           is_candidate_modulator,NEPC_viper_score,is_a_known_marker,
                                           is_c_mrs,gene_human,cindy_pvalue_integrated,
                                           cindy_neg_log_pvalue_integrated,
                                           modulator_ranking,dge_FDR=FDR,dge_logFC=logFC)
  
  my_tibble <- my_tibble %>% arrange(cis_fisher_pvalue,is_cis_gene,is_diff_expr,modulator_ranking)
  
  my_tibble$F <- NULL
  my_tibble$dge_logFC <- round(my_tibble$dge_logFC,3)
  my_tibble$NEPC_viper_score <- round(my_tibble$NEPC_viper_score,3)
  return(my_tibble)
}
saveNEPCSignaturesTables <- function(my_tibble){
  print_msg_info(">>> >> Saving Table to Excel")
  writexl::write_xlsx(my_tibble, file.path("experiments/integrative-analysis/reports/nepc-signatures-table-full.xlsx") )
  writexl::write_xlsx(my_tibble %>% slice_head(n = 100) , file.path("experiments/integrative-analysis/reports/nepc-signatures-table-reduced.xlsx") )
  writexl::write_xlsx(my_tibble %>% filter(is_cis_gene=="Yes") , file.path("experiments/integrative-analysis/reports/nepc-signatures-table-all-cis.xlsx") )
  writexl::write_xlsx(x <- my_tibble %>% filter(is_cis_gene=="Yes",is_diff_expr=="Yes",is_candidate_modulator=="Yes") , file.path("experiments/integrative-analysis/reports/nepc-signatures-table-all-3-conditions.xlsx") )
}
# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
# -------------------------------------------------------------------------------------------
# ----------------------------------- ADD REGULATORS FUNC -----------------------------------
# -------------------------------------------------------------------------------------------
# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
addRegulatorTypeColumnToMyTibble <- function(my_tibble, tfs, cotfs, sig){
  my_tibble$regulator_type <- NA
  my_tibble$regulator_type <- ifelse( my_tibble$gene_human %in% tfs$X1 , "TF" , my_tibble$regulator_type )
  my_tibble$regulator_type <- ifelse( my_tibble$gene_human %in% cotfs$X1 , "coTF" , my_tibble$regulator_type )
  my_tibble$regulator_type <- ifelse( my_tibble$gene_human %in% sig$X1 , "SIG" , my_tibble$regulator_type )
  return(my_tibble)
}
# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
# -------------------------------------------------------------------------------------------
# ------------------------------ SAVE INTEGRATIVE ANALYSIS FUNC -----------------------------
# -------------------------------------------------------------------------------------------
# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
saveSbCisNepcViperCindyIntegrativeAnalysisToExcel <- function(my_tibble){
  filename <- "experiments/integrative-analysis/processed_data/tibble-of-sb-cis-nepc-viper-cindy-integrative-analysis.rds"
  saveRDS(my_tibble,filename)
}
# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
# -------------------------------------------------------------------------------------------
# ---------------------------------- FISCHER INTEGRATION FUNCS ----------------------------------
# -------------------------------------------------------------------------------------------
# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&

# ------------------------------------- HELPER FUNCTIONS ------------------------------------
my_p_integration_func <- function( my_pvalues ) {
  pchisq(q = (-2*sum(log(my_pvalues))), df = (2*length(my_pvalues)), lower.tail = FALSE)	
}

# -------------------------------------- MAIN FUNCTIONS -------------------------------------
performFisherIntegration <- function(my_tibble){
  print_msg_info(">>> Fisher's Integration")
  my_tibble$integration_cis_p <- my_tibble$cis_fisher_pvalue
  my_tibble$integration_viper_p <- pnorm( abs( my_tibble$NEPC_viper_score ) , lower.tail = FALSE )
  my_tibble$integration_dge_p <- my_tibble$dge_FDR
  my_tibble$integration_cindy_p <- my_tibble$cindy_pvalue_integrated
  
  my_tibble$integration_cis_p <- ifelse( is.na(my_tibble$integration_cis_p) | is.nan(my_tibble$integration_cis_p) , 1 , my_tibble$integration_cis_p )
  my_tibble$integration_viper_p <- ifelse( is.na(my_tibble$integration_viper_p) | is.nan(my_tibble$integration_viper_p) , 1 , my_tibble$integration_viper_p )
  my_tibble$integration_dge_p <- ifelse( is.na(my_tibble$integration_dge_p) | is.nan(my_tibble$integration_dge_p) , 1 , my_tibble$integration_dge_p )
  my_tibble$integration_cindy_p <- ifelse( is.na(my_tibble$integration_cindy_p) | is.nan(my_tibble$integration_cindy_p) , 1 , my_tibble$integration_cindy_p )
  my_tibble$integration_cindy_neg_log_p <- ifelse( is.na(my_tibble$integration_cindy_p) | is.nan(my_tibble$integration_cindy_p) , 1 , my_tibble$cindy_neg_log_pvalue_integrated )
  
  my_mat <- my_tibble %>% 
    dplyr::select(integration_cis_p,integration_viper_p,integration_dge_p) %>%
    as.matrix()
  
  # print_msg_warn("*** Missing Integration with CINDy score ***")	
  
  my_pvalues <- apply( my_mat , 1 , my_p_integration_func )
  
  my_tibble$final_integration_pvalues <- my_pvalues
  
  my_tibble <- my_tibble %>% dplyr::select(gene,final_integration_pvalues,everything()) %>% arrange(final_integration_pvalues)
  return(my_tibble)
}
saveFisherIntegrationTableToExcel <- function(my_tibble){
  my_tibble$cindy_pvalue_integrated <- NULL
  my_tibble$cindy_neg_log_pvalue_integrated <- NULL
  writexl::write_xlsx(my_tibble, file.path("experiments/integrative-analysis/reports/nepc-signatures-table-final-fisher-integration-table.xlsx") )
  write_tsv(my_tibble, file.path("experiments/integrative-analysis/reports/nepc-signatures-table-final-fisher-integration-table.tsv") )
  write_tsv(my_tibble, file.path("experiments/integrative-analysis/processed_data/nepc-signatures-table-final-fisher-integration-table.tsv") )
  
  writexl::write_xlsx(my_tibble %>% filter(is_cis_gene=="Yes",is_diff_expr=="Yes",is_candidate_modulator=="Yes")  , 
                      file.path("experiments/integrative-analysis/reports/top-hits-for-validation-final-fisher-integration-table.xlsx") )
  writexl::write_xlsx(my_tibble %>% filter(is_cis_gene=="Yes",is_diff_expr=="Yes",is_candidate_modulator=="Yes")  , 
                      file.path("experiments/integrative-analysis/processed_data/top-hits-for-validation-final-fisher-integration-table.xlsx") )	
}

# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
# -------------------------------------------------------------------------------------------
# ----------------------------------- PLOTTING CCLE FUNCS -----------------------------------
# -------------------------------------------------------------------------------------------
# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
## Plotting  log2(TPM)
plot_log2_tpm <- function(my_tibble, ccle_prad.log2tpm){
  my_tibble.filtered <- my_tibble %>% 
    filter(is_cis_gene=="Yes",is_diff_expr=="Yes",is_candidate_modulator=="Yes") %>% 
    arrange(final_integration_pvalues)
  sorted_pca_cline.names <- c("NCIH660_PROSTATE","LNCAPCLONEFGC_PROSTATE","DU145_PROSTATE","22RV1_PROSTATE",
                              "MDAPCA2B_PROSTATE","PC3_PROSTATE","PRECLH_PROSTATE","VCAP_PROSTATE")
  index <- my_tibble$is_cis_gene == "Yes" & my_tibble$is_diff_expr == "Yes" & my_tibble$is_candidate_modulator == "Yes"
  index <- rownames(ccle_prad.log2tpm) %in% my_tibble$gene_human[ index ]
  mat <- ccle_prad.log2tpm[ index , ]
  mat <- mat[ match( my_tibble.filtered$gene_human , rownames(mat) ) , ]
  
  mat <- mat[ , match( sorted_pca_cline.names , colnames(mat) ) ]
  colnames(mat) <- gsub( "(.*)_(.*)" , "\\1" , colnames(mat) )
  
  cis2validate_on_cell_lines_hm <- Heatmap( mat , 
                                            col = viridis(10) ,
                                            cluster_rows = FALSE , cluster_columns = FALSE )																							
  
  p <- draw(cis2validate_on_cell_lines_hm) # , heatmap_legend_side = "bottom" )
  
  pdf( file.path( "experiments/integrative-analysis/reports/log2tpm-of-cis2validate-on-prostate-cancer-cell-lines.pdf")  ,
       width = 4 , height = 8 )
  print(p)
  dev.off() 
}
## Plotting scaled log2(TPM)
plot_scaled_log2_tpm <- function(my_tibble, ccle_prad.log2tpm){
  sorted_pca_cline.names <- c("NCIH660_PROSTATE","LNCAPCLONEFGC_PROSTATE","DU145_PROSTATE","22RV1_PROSTATE",
                              "MDAPCA2B_PROSTATE","PC3_PROSTATE","PRECLH_PROSTATE","VCAP_PROSTATE")
  my_tibble.filtered <- my_tibble %>% 
    filter(is_cis_gene=="Yes",is_diff_expr=="Yes",is_candidate_modulator=="Yes") %>% 
    arrange(final_integration_pvalues)
  index <- my_tibble$is_cis_gene == "Yes" & my_tibble$is_diff_expr == "Yes" & my_tibble$is_candidate_modulator == "Yes"
  index <- rownames(ccle_prad.log2tpm) %in% my_tibble$gene_human[ index ]
  mat <- ccle_prad.log2tpm[ index , ]
  scaled_mat <- t(scale(t(mat)))
  scaled_mat <- scaled_mat[ , match( sorted_pca_cline.names , colnames(scaled_mat) ) ]
  scaled_mat <- scaled_mat[ match( my_tibble.filtered$gene_human , rownames(scaled_mat) ) , ]
  colnames(scaled_mat) <- gsub( "(.*)_(.*)" , "\\1" , colnames(scaled_mat) )
  col_fun = colorRamp2(c(-2, 0, +2), c("blue", "white", "red"))
  cis2validate_on_cell_lines_hm <- Heatmap( scaled_mat , 
                                            col = col_fun ,
                                            cluster_rows = FALSE , cluster_columns = FALSE )
  
  p <- draw(cis2validate_on_cell_lines_hm) # , heatmap_legend_side = "bottom" )
  
  pdf( file.path( "experiments/integrative-analysis/reports/scaled-log2tpm-of-cis2validate-on-prostate-cancer-cell-lines.pdf")  , width = 4 , height = 8 )
  print(p)
  dev.off()
}
## Plotting scaled log2(TPM)
plot_scaled_log2_tpm_all_ccle <- function(my_tibble, ccle_prad.log2tpm, ccle.log2tpm){
  sorted_pca_cline.names <- c("NCIH660_PROSTATE","LNCAPCLONEFGC_PROSTATE","DU145_PROSTATE","22RV1_PROSTATE",
                              "MDAPCA2B_PROSTATE","PC3_PROSTATE","PRECLH_PROSTATE","VCAP_PROSTATE")
  index <- my_tibble$is_cis_gene == "Yes" & my_tibble$is_diff_expr == "Yes" & my_tibble$is_candidate_modulator == "Yes"
  index <- rownames(ccle_prad.log2tpm) %in% my_tibble$gene_human[ index ]
  
  my_tibble <- my_tibble %>% arrange(final_integration_pvalues)
  index <- rownames(ccle_prad.log2tpm) %in% my_tibble$gene_human[ index ]
  mat <- ccle.log2tpm[ index , ]
  scaled_mat <- t(scale(t(mat)))
  scaled_mat <- scaled_mat[ , match( sorted_pca_cline.names , colnames(scaled_mat) ) ]
  scaled_mat <- scaled_mat[ match( my_tibble$gene_human , rownames(scaled_mat) ) , ]
  colnames(scaled_mat) <- gsub( "(.*)_(.*)" , "\\1" , colnames(scaled_mat) )
  col_fun = colorRamp2(c(-2, 0, +2), c("blue", "white", "red"))		
  cis2validate_on_cell_lines_hm <- Heatmap( scaled_mat , 
                                            col = col_fun ,
                                            cluster_rows = FALSE , cluster_columns = FALSE )
  
  p <- draw(cis2validate_on_cell_lines_hm) # , heatmap_legend_side = "bottom" )
  
  pdf( file.path( "experiments/integrative-analysis/reports/scaled-log2tpm-of-cis2validate-on-cancer-cell-lines.pdf")  , width = 4 , height = 8 )
  print(p)
  dev.off()  				
}
## Plotting scaled log2(TPM) for all cis (not only the ones selected for validation)
plot_scaled_log2_tpm_all_cis <- function(my_tibble, ccle_prad.log2tpm){
  sorted_pca_cline.names <- c("NCIH660_PROSTATE","LNCAPCLONEFGC_PROSTATE","DU145_PROSTATE","22RV1_PROSTATE",
                              "MDAPCA2B_PROSTATE","PC3_PROSTATE","PRECLH_PROSTATE","VCAP_PROSTATE")
  my_tibble.filtered <- my_tibble %>% 
    filter(is_cis_gene=="Yes") %>% 
    arrange(cis_fisher_pvalue)
  index <- rownames(ccle_prad.log2tpm) %in% my_tibble.filtered$gene_human
  mat <- ccle_prad.log2tpm[ index , ]
  scaled_mat <- t(scale(t(mat)))
  scaled_mat <- scaled_mat[ , match( sorted_pca_cline.names , colnames(scaled_mat) ) ]
  scaled_mat <- scaled_mat[ match( my_tibble.filtered$gene_human , rownames(scaled_mat) ) , ]
  colnames(scaled_mat) <- gsub( "(.*)_(.*)" , "\\1" , colnames(scaled_mat) )
  col_fun = colorRamp2(c(-2, 0, +2), c("blue", "white", "red"))		
  cis2validate_on_cell_lines_hm <- Heatmap( scaled_mat , 
                                            col = col_fun ,
                                            cluster_rows = FALSE , cluster_columns = FALSE )
  
  p <- draw(cis2validate_on_cell_lines_hm) # , heatmap_legend_side = "bottom" )
  
  pdf( file.path( "experiments/integrative-analysis/reports/scaled-log2tpm-of-allCis-on-prostate-cancer-cell-lines.pdf")  , width = 3.5 , height = 42 )
  print(p)
  dev.off()
}
## Plotting scaled log2(TPM) for some markers
plot_scaled_log2_tpm_markers <- function(my_tibble, ccle_prad.log2tpm){
  sorted_pca_cline.names <- c("NCIH660_PROSTATE","LNCAPCLONEFGC_PROSTATE","DU145_PROSTATE","22RV1_PROSTATE",
                              "MDAPCA2B_PROSTATE","PC3_PROSTATE","PRECLH_PROSTATE","VCAP_PROSTATE")
  index <- rownames(ccle_prad.log2tpm) %in% c("KLK3","FOXM1","CENPF","CHGA","SYP","AR","NR3C1","MYC","MYCN","EZH2")
  mat <- ccle_prad.log2tpm[ index , ]
  scaled_mat <- t(scale(t(mat)))
  scaled_mat <- scaled_mat[ , match( sorted_pca_cline.names , colnames(scaled_mat) ) ]
  colnames(scaled_mat) <- gsub( "(.*)_(.*)" , "\\1" , colnames(scaled_mat) )
  col_fun = colorRamp2(c(-2, 0, +2), c("blue", "white", "red"))		
  cis2validate_on_cell_lines_hm <- Heatmap( scaled_mat , 
                                            col = col_fun ,
                                            cluster_rows = FALSE , cluster_columns = FALSE )
  
  p <- draw(cis2validate_on_cell_lines_hm) # , heatmap_legend_side = "bottom" )
  
  pdf( file.path( "experiments/integrative-analysis/reports/scaled-log2tpm-of-markers-on-prostate-cancer-cell-lines.pdf")  , width = 4 , height = 6 )
  print(p)
  dev.off() 		
}
## Plotting scaled log2(TPM) for some markers 
plot_scaled_log2_tpm_markers_on_all_ccle <- function(my_tibble, ccle_prad.log2tpm, ccle.log2tpm){
  sorted_pca_cline.names <- c("NCIH660_PROSTATE","LNCAPCLONEFGC_PROSTATE","DU145_PROSTATE","22RV1_PROSTATE",
                              "MDAPCA2B_PROSTATE","PC3_PROSTATE","PRECLH_PROSTATE","VCAP_PROSTATE")
  index <- rownames(ccle_prad.log2tpm) %in% c("KLK3","FOXM1","CENPF","CHGA","SYP","AR","NR3C1","MYC","MYCN","EZH2")
  mat <- ccle.log2tpm[ index , ]
  scaled_mat <- t(scale(t(mat)))
  scaled_mat <- scaled_mat[ , match( sorted_pca_cline.names , colnames(scaled_mat) ) ]
  colnames(scaled_mat) <- gsub( "(.*)_(.*)" , "\\1" , colnames(scaled_mat) )
  col_fun = colorRamp2(c(-2, 0, +2), c("blue", "white", "red"))		
  cis2validate_on_cell_lines_hm <- Heatmap( scaled_mat , 
                                            col = col_fun ,
                                            cluster_rows = FALSE , cluster_columns = FALSE )
  
  p <- draw(cis2validate_on_cell_lines_hm) # , heatmap_legend_side = "bottom" )
  
  pdf( file.path( "experiments/integrative-analysis/reports/scaled-log2tpm-of-markers-on-cancer-cell-lines.pdf")  , width = 4 , height = 6 )
  print(p)
  dev.off()
}

## Integrative Analysis ----
print_msg_info(">>> Integrative Analysis")

## > Loading CINDy Table ----
print_msg_info(">>> >> Loading CINDy Table")
cindy_table <- readRDS(filepath_cindy_table_su2c_plus_tcga_merged_rds)
cindy_table <- preprocessCindyTable(cindy_table)

## > Loading Gene Expression and VIPER analysis ----
print_msg_info(">>> >> Loading Gene Expression and VIPER analysis")
# my_tibble <- readRDS("data/mouse-analysis/NPp53_NEPC_effect_on_SB_activated-tibble-manual-genes-curation.rds")
my_tibble <- readRDS(filepath_NPp53_NEPC_effect_on_SB_activated_tibble_rds)
# my_tibble <- readRDS("experiments/all-samples-sb-rna-seq-analysis/reports/edgeR/NPp53_SB_effect_on_NEPC/NPp53_SB_effect_on_NEPC-tibble.rds")
vpmat <- readRDS(filepath_vpmat_rds)
# vpmat_old <- readRDS("data_old/mouse-analysis/vpmat.rds")
my_tibble <- preprocessMyTibble(my_tibble, vpmat)

## > Gathering CIS-assiciated genes ----
print_msg_info(">>> >> Gathering CIS-assiciated genes")
print_msg_warn(">>> *** This List is not derived from NE-only samples ***")
cis_table <- read_csv(filepath_cis_all_nr_prostate_0_001_csv, skip = 1)
cis_table <- preprocessCisTable(cis_table)
## > Gathering SB RNA-Seq Metadata ----
print_msg_info(">>> >> > Gathering SB RNA-Seq Metadata")
sb_excel <- readxl::read_xlsx(filepath_SB_List_CIS_and_RNAseq_2020_OCT_20_xlsx, sheet = 2)
sb_excel <- preprocessSbExcel(sb_excel)
## Separating rows by making a row for each Library and copying the rest of the columns
saveCisTablePerLibrary(cis_table, sb_excel)

## > Integrating DATA ----
my_tibble <- integrateCindyTableAndCisTableIntoMyTibble(my_tibble, cindy_table, cis_table)

## > Saving Table to Excel ----
tfs <- read_csv(filepath_tfs_csv,col_names = F)
cotfs <- read_csv(filepath_cotfs_csv,col_names = F)
sig <- read_csv(filepath_sig_csv,col_names = F)
my_tibble <- reformatMyTibbleToSaveNEPCSignatures(my_tibble)
saveNEPCSignaturesTables(my_tibble)
my_tibble <- addRegulatorTypeColumnToMyTibble(my_tibble, tfs, cotfs, sig)
saveSbCisNepcViperCindyIntegrativeAnalysisToExcel(my_tibble)
my_tibble <- performFisherIntegration(my_tibble)
saveFisherIntegrationTableToExcel(my_tibble)

## > CIS on CCLE  ----
print_msg_info(">>> >> CIS on CCLE")
ccle.tpm <- readRDS(filepath_ccle_tpm_symbols_rds)
ccle_annotations <- read_delim( file = filepath_Cell_lines_annotations_20181226_txt, delim = "\t" )
ccle_prad.log2tpm <- preprocessCCLE(ccle.tpm, ccle_annotations)
plot_log2_tpm(my_tibble, ccle_prad.log2tpm)
plot_scaled_log2_tpm(my_tibble, ccle_prad.log2tpm)
plot_scaled_log2_tpm_all_ccle(my_tibble, ccle_prad.log2tpm, ccle.log2tpm = log2(ccle.tpm+1))
plot_scaled_log2_tpm_all_cis(my_tibble, ccle_prad.log2tpm)
plot_scaled_log2_tpm_markers(my_tibble, ccle_prad.log2tpm)
plot_scaled_log2_tpm_markers_on_all_ccle(my_tibble, ccle_prad.log2tpm, ccle.log2tpm = log2(ccle.tpm+1))
