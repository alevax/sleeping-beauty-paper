suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(ggsurvfit)))
suppressMessages(suppressWarnings(library(ComplexHeatmap)))

# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
# -----------------------------------------------------------------------------------
# ------------------------------- LOAD DATA FUNCTIONS -------------------------------
# -----------------------------------------------------------------------------------
# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
GESTransform <- function(dat.mat, na.rm = TRUE){
  # generate GES
  ges.mat <- t(apply(dat.mat, 1, function(x) {
    (x - mean(x, na.rm = na.rm)) / sd(x, na.rm = na.rm)
  }))
  ges.mat
}
getCisTopHitsAndVPScores <- function(my_top_hits_table){
  cis_top_hits <- my_top_hits_table$gene_human
  cis_top_hits_vpscores <- my_top_hits_table$NEPC_viper_score
  names(cis_top_hits_vpscores) <- cis_top_hits
  cis_top_hits_vpscores <- cis_top_hits_vpscores[!is.na(cis_top_hits_vpscores)]
  cis_top_hits <- names(cis_top_hits_vpscores)
  my_cis_top_hits_list <- list(cis_top_hits, cis_top_hits_vpscores)
  names(my_cis_top_hits_list) <- c("cis_top_hits", "cis_top_hits_vpscores")
  return(my_cis_top_hits_list)
}
loadSu2CECMetadata <- function(){
  SU2CEC_metadata <- readRDS("experiments/survival-analysis/processed_data/metadata/SU2C_EastCoast_metadata_forSA.rds") %>%
    dplyr::select("Sample_Identifier", "Neuroendocrine_Features", "Overall_Survival_Status", "Overall_Survival_Days") %>%
    dplyr::filter(!is.na(Overall_Survival_Days)) %>%
    dplyr::mutate(time = Overall_Survival_Days) %>%
    dplyr::mutate(status = (Overall_Survival_Status == "1:DECEASED")) %>%
    dplyr::filter(!is.na(Neuroendocrine_Features))
  return(SU2CEC_metadata)
}
loadSu2CWCMetadata <- function(){
  SU2CWC_metadata <- readRDS("experiments/oncomatch-analysis/processed_data/metadata/SU2C_WestCoast_metadata_forOMPlot.rds")
  JCO_2018 <- as_tibble(read.csv("data/SU2C_WestCoast/WCDT/clinical_metadata/2021_11_21_JCO_2018.txt", sep = "\t"))
  # Update sample IDs to keep them consistent
  JCO_2018$sample_id[JCO_2018$sample_id == "DTB-098-PRO"] <- "DTB-098-PRO2"
  JCO_2018$sample_id[JCO_2018$sample_id == "DTB-092-BL"] <- "DTB-092-BL-2018"
  JCO_2018$sample_id[JCO_2018$sample_id == "DTB-055-PRO2"] <- "DTB-055-PRO"
  JCO_2018 <- JCO_2018 %>%
    dplyr::select("sample_id", "OS_duration", "OS_event_.dead.") %>%
    dplyr::rename("Sample_Identifier" = "sample_id",
                  "Overall_Survival_Status" = "OS_event_.dead.",
                  "Overall_Survival_Days" = "OS_duration") %>%
    dplyr::mutate("time" = as.numeric(Overall_Survival_Days)) %>%
    dplyr::mutate("status" = as.logical(Overall_Survival_Status))
  SU2CWC_metadata <- SU2CWC_metadata %>%
    dplyr::select(Sample_Identifier, Neuroendocrine_Features) %>%
    left_join(JCO_2018, by = c("Sample_Identifier" = "Sample_Identifier")) %>%
    dplyr::mutate(status = (Overall_Survival_Status == "Death")) %>%
    dplyr::filter(!is.na(Overall_Survival_Days)) %>%
    dplyr::filter(!is.na(Neuroendocrine_Features))
  return(SU2CWC_metadata)
}
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

# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
# -----------------------------------------------------------------------------------
# ------------------------ INDIVIDUAL SURVIVAL ANALYSIS FUNCS -----------------------
# -----------------------------------------------------------------------------------
# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
individual_survival_by_gene_activity <- function(SU2CEC_metadata, gene_i, vpscores, out_dir, thresh_method = "quantile"){
  # Filter the data down to significant activity
  if(thresh_method == "quantile"){
    SU2CEC_metadata_2 <- SU2CEC_metadata %>%
      dplyr::filter(get(gene_i) > quantile(SU2CEC_metadata[[gene_i]], 0.75) | get(gene_i) < quantile(SU2CEC_metadata[[gene_i]], 0.25))
  } else {
    SU2CEC_metadata_2 <- SU2CEC_metadata %>%
      dplyr::filter(get(gene_i) >= 1.96 | get(gene_i) <= -1.96)
  }
  
  if(nrow(SU2CEC_metadata_2) < 2){
    return(NA)
  }

  # Calculate the p-value
  gene_i_sign <- ifelse(sign(vpscores[gene_i]) == 1, "+", "-")
  gene_i_antisign <- ifelse(sign(vpscores[gene_i]) == 1, "-", "+")
  SU2CEC_metadata_2[[gene_i]] <- as.integer(SU2CEC_metadata_2[[gene_i]] > 0)
  gene_i_coxph <- survival::coxph(Surv(time, status) ~ get(gene_i), data = SU2CEC_metadata_2)
  p.value <- gtsummary::tbl_regression(gene_i_coxph, exp = TRUE)$table_body$p.value
  
  if(is.na(p.value)){
    return(NA)
  }
  if(p.value <= 0.05){
    out_dir <- paste0(out_dir, "stats/sig/")
  } else {
    out_dir <- paste0(out_dir, "stats/nonSig/")
  }
  
  # Save the txt
  my_filename <- paste0(out_dir, gene_i, "_stats.txt")
  write(capture.output(print(gene_i_coxph)), my_filename)
  
  # Save the survival plot
  SU2CEC_metadata_2[[gene_i]] <- as.factor(SU2CEC_metadata_2[[gene_i]])
  levels(SU2CEC_metadata_2[[gene_i]]) <- c(
    paste0(gene_i, gene_i_antisign),
    paste0(gene_i, gene_i_sign)
  )
  
  survival_plot <- survfit2(Surv(time, status) ~ get(gene_i), data = SU2CEC_metadata_2) %>% 
    ggsurvfit() +
    labs(
      x = "Days",
      y = "Overall survival probability"
    ) + 
    add_confidence_interval() +
    scale_color_manual(values = c("blue", "red")) +
    scale_fill_manual(values = c("blue", "red")) +
    add_risktable() +
    add_pvalue(location = "annotation", size = 5)
  my_res_list <- list(survival_plot, p.value)
  names(my_res_list) <- c("survival_plot", "p.value")
  return(my_res_list)
}

getSurvivalForGeneList <- function(SU2CEC_metadata,
                                   SU2CEC_vpmat,
                                   gene_list,
                                   vpscores,
                                   out_dir,
                                   out_name,
                                   thresh_method = "quantile"){
  dir.create(paste0(out_dir, "stats/"), showWarnings = FALSE)
  dir.create(paste0(out_dir, "stats/sig/"), showWarnings = FALSE)
  dir.create(paste0(out_dir, "stats/nonSig/"), showWarnings = FALSE)
  res_df <- data.frame(matrix(NA, nrow = length(gene_list), ncol = 2))
  colnames(res_df) <- c("gene", "p_value")
  rownames(res_df) <- gene_list
  res_df$gene <- gene_list
  survival_plot_list_signif <- list()
  survival_plot_list_nonSignif <- list()
  for(i in 1:length(gene_list)){
    gene_i <- gene_list[i]
    if(gene_i %in% rownames(SU2CEC_vpmat)){
      SU2CEC_metadata <- cbind(SU2CEC_metadata, SU2CEC_vpmat[gene_i, SU2CEC_metadata$Sample_Identifier]*sign(vpscores[gene_i]))
      colnames(SU2CEC_metadata)[ncol(SU2CEC_metadata)] <- gene_i
      my_res_list <- individual_survival_by_gene_activity(SU2CEC_metadata, gene_i, vpscores, out_dir, thresh_method)
      if(!any(is.na(my_res_list))){
        res_df[gene_i, "p_value"] <- my_res_list[["p.value"]]
        gene_i_plot_list <- list(my_res_list[["survival_plot"]])
        names(gene_i_plot_list) <- gene_i
        if(my_res_list[["p.value"]] <= 0.05){
          survival_plot_list_signif <- c(survival_plot_list_signif, gene_i_plot_list)
        } else {
          survival_plot_list_nonSignif <- c(survival_plot_list_nonSignif, gene_i_plot_list)
        }
      }
    }
  }
  # Save CSVs of p-values
  saveRDS(object = res_df, file = paste0(out_dir, out_name, "_pvals_df.rds"))
  res_df$p_value <- round(res_df$p_value, 4)
  res_df <- res_df %>%
    as_tibble() %>%
    arrange(p_value)
  write.csv(x = res_df, file = paste0(out_dir, out_name, "_pvals_df.csv"))
  res_df_signif <- res_df %>%
    filter(p_value <= 0.05)
  write.csv(x = res_df_signif, file = paste0(out_dir, out_name, "_pvals_df_signif.csv"))
  # Save PDFs of survival plots
  my_filename <- paste0(out_dir, out_name, "_survival_signif.pdf")
  pdf(file = my_filename,
      onefile = TRUE,
      width = 10,
      height = 12)
  print(survival_plot_list_signif[res_df_signif$gene])
  dev.off()
  my_filename <- paste0(out_dir, out_name, "_survival_nonSignif.pdf")
  pdf(file = my_filename,
      onefile = TRUE,
      width = 10,
      height = 12)
  print(survival_plot_list_nonSignif[setdiff(res_df$gene, res_df_signif$gene)])
  dev.off()
}

# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
# -----------------------------------------------------------------------------------
# ------------------------- PAIRWISE SURVIVAL ANLAYSIS FUNCS ------------------------
# -----------------------------------------------------------------------------------
# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
pairwise_survival_by_gene_activity <- function(SU2CEC_metadata_2,
                                               gene_i,
                                               gene_j,
                                               vpscores,
                                               out_dir,
                                               genes_type = "cis_top_hits",
                                               thresh_method = "quantile"){
  # Filter the data down to significant activity
  if(thresh_method == "quantile"){
    SU2CEC_metadata_3 <- SU2CEC_metadata_2 %>%
      dplyr::mutate("genes_i_j_active" = (get(gene_i) >= quantile(SU2CEC_metadata_2[[gene_i]], 0.75) & get(gene_j) >= quantile(SU2CEC_metadata_2[[gene_j]], 0.75))) %>%
      dplyr::mutate("genes_active" = as.integer(genes_i_j_active))
  } else {
    SU2CEC_metadata_3 <- SU2CEC_metadata_2 %>%
      dplyr::mutate("genes_i_j_active" = (get(gene_i) >= 1.96 & get(gene_j) <= -1.96)) %>%
      dplyr::mutate("genes_active" = as.integer(genes_i_j_active))
  }
  
  
  if(sum(SU2CEC_metadata_3$genes_i_j_active) < 2){
    return(NA)
  }
  
  genes_names_for_files <- paste0(gene_i, "&", gene_j)
  gene_i_sign <- ifelse(sign(vpscores[gene_i]) == 1, "+", "-")
  gene_j_sign <- ifelse(sign(vpscores[gene_j]) == 1, "+", "-")
  SU2CEC_metadata_3[["genes_active"]] <- as.factor(SU2CEC_metadata_3[["genes_active"]])
  levels(SU2CEC_metadata_3[["genes_active"]]) <- c(
    paste0(gene_i, "+/-", " & ", gene_j, "-/+"),
    paste0(gene_i, gene_i_sign, " & ", gene_j, gene_j_sign)
  )
  
  # Calculate the p-value
  gene_i_coxph <- survival::coxph(Surv(time, status) ~ genes_active, data = SU2CEC_metadata_3)
  gene_i_coxph_plot <- gtsummary::tbl_regression(gene_i_coxph, exp = TRUE) 
  p.value <- gene_i_coxph_plot$table_body$p.value[3]
  if(is.na(p.value)){
    return(NA)
  }
  if(p.value <= 0.05){
    out_dir <- paste0(out_dir, "stats/sig/")
  } else {
    out_dir <- paste0(out_dir, "stats/nonSig/")
  }
  
  # Save the txt
  my_filename <- paste0(out_dir, genes_names_for_files, "_stats.txt")
  write(capture.output(print(gene_i_coxph)), my_filename)
  # Save plot as PNG
  if(genes_type == "cis_top_hits" || p.value <= 0.05){
    survival_plot <- survfit2(Surv(time, status) ~ genes_active, data = SU2CEC_metadata_3) %>%
      ggsurvfit() +
      labs(
        x = "Days",
        y = "Overall survival probability"
      ) +
      add_confidence_interval() +
      scale_color_manual(values = c("blue", "red")) + 
      scale_fill_manual(values = c("blue", "red")) +
      add_risktable() +
      add_pvalue(location = "annotation", size = 5)
    my_res_list <- list(survival_plot, p.value)
    names(my_res_list) <- c("survival_plot", "p.value")
  } else {
    my_res_list <- list(p.value)
    names(my_res_list) <- c("p.value")
  }
  return(my_res_list)
}
performPairChecksAndSaveResults <- function(SU2CEC_metadata,
                                            SU2CEC_vpmat,
                                            genes_names,
                                            vpscores,
                                            out_dir,
                                            save_plots = TRUE,
                                            genes_type = "cis_top_hits",
                                            thresh_method = "quantile"){
  n_top_hits <- length(genes_names)
  
  dir.create(out_dir, showWarnings = FALSE)
  dir.create(paste0(out_dir, "stats/"), showWarnings = FALSE)
  dir.create(paste0(out_dir, "stats/sig/"), showWarnings = FALSE)
  dir.create(paste0(out_dir, "stats/nonSig/"), showWarnings = FALSE)
  
  res_mat <- matrix(NA, nrow = n_top_hits, ncol = n_top_hits)
  colnames(res_mat) <- rownames(res_mat) <- genes_names
  pb = txtProgressBar(min = 0, max = n_top_hits, initial = 0, style = 3)
  
  res_df <- data.frame(matrix(NA, nrow = 0, ncol = 3))
  colnames(res_df) <- c("gene_1", "gene_2", "p_value")
  # Make sure the signs are correct in the plots
  
  n_samples <- length(SU2CEC_metadata$Sample_Identifier)
  
  signif_plot_list <- list()
  if(genes_type == "cis_top_hits") nonSignif_plot_list <- list()
  
  for(i in 1:n_top_hits){
    for(j in 1:n_top_hits){
      if(i > j){
        gene_i <- genes_names[i]
        gene_j <- genes_names[j]
        SU2CEC_metadata_2 <- cbind(SU2CEC_metadata,
                                   SU2CEC_vpmat[gene_i, SU2CEC_metadata$Sample_Identifier]*sign(vpscores[gene_i]))
        colnames(SU2CEC_metadata_2)[ncol(SU2CEC_metadata_2)] <- gene_i
        SU2CEC_metadata_2 <- cbind(SU2CEC_metadata_2,
                                   SU2CEC_vpmat[gene_j, SU2CEC_metadata$Sample_Identifier]*sign(vpscores[gene_j]))
        colnames(SU2CEC_metadata_2)[ncol(SU2CEC_metadata_2)] <- gene_j
        my_res_list <- pairwise_survival_by_gene_activity(SU2CEC_metadata_2,
                                                          gene_i,
                                                          gene_j,
                                                          vpscores = vpscores,
                                                          out_dir,
                                                          genes_type,
                                                          thresh_method)
        if(!any(is.na(my_res_list))){
          pval <- my_res_list[["p.value"]]
          res_mat[gene_i, gene_j] <- pval
          res_mat[gene_j, gene_i] <- pval
          
          genes_i_j_plot_list <- list(my_res_list[["survival_plot"]])
          names(genes_i_j_plot_list) <- paste0(gene_i, "&", gene_j)
          if(pval <= 0.05){
            signif_plot_list <- c(signif_plot_list, genes_i_j_plot_list)
          } else if(genes_type == "cis_top_hits"){
            nonSignif_plot_list <- c(nonSignif_plot_list, genes_i_j_plot_list)
          }
        } else {
          pval <- NA
        }
        res_df <- rbind(res_df, matrix(c(gene_i, gene_j, pval), nrow = 1, ncol = 3))
        setTxtProgressBar(pb, i)
      }
    }
  }
  close(pb)
  
  # Use the script above to return p-values and put these into a matrix, which you can then turn into a heatmap
  res_mat_binary <- res_mat <= 0.05
  res_mat_binary <- res_mat_binary[names(sort(rowSums(res_mat_binary, na.rm=TRUE), decreasing = TRUE)),
                                   names(sort(colSums(res_mat_binary, na.rm=TRUE), decreasing = TRUE))]
  p <- draw(pheatmap(res_mat_binary,
                     cluster_rows = FALSE,
                     cluster_cols = FALSE,
                     legend = FALSE))
  
  if(length(genes_names) > 50){
    pdf(file = paste0(out_dir, "pair_matches_summary.pdf"),
        onefile = FALSE,
        width = 18,
        height = 15)
  } else {
    pdf(file = paste0(out_dir, "pair_matches_summary.pdf"),
        onefile = FALSE,
        width = 12,
        height = 10)
  }
  print(p)
  dev.off()
  
  p_bw <- draw(pheatmap(res_mat_binary,
                        cluster_rows = FALSE,
                        cluster_cols = FALSE,
                        color = c("white", "black"),
                        legend = FALSE))
  if(length(genes_names) > 50){
    pdf(file = paste0(out_dir, "pair_matches_summary_bw.pdf"),
        onefile = FALSE,
        width = 18,
        height = 15)
  } else {
    pdf(file = paste0(out_dir, "pair_matches_summary_bw.pdf"),
        onefile = FALSE,
        width = 12,
        height = 10)
  }
  print(p_bw)
  dev.off()
  
  res_mat_logp <- log(res_mat, base = 0.05)
  res_mat_logp <- res_mat_logp[names(sort(rowSums(res_mat_logp, na.rm=TRUE), decreasing = TRUE)),
                               names(sort(colSums(res_mat_logp, na.rm=TRUE), decreasing = TRUE))]
  res_mat_logp <- res_mat_logp[names(sort(rowSums(res_mat_binary, na.rm=TRUE), decreasing = TRUE)),
                               names(sort(colSums(res_mat_binary, na.rm=TRUE), decreasing = TRUE))]
  
  p2 <- draw(pheatmap(res_mat_logp,
                      # color = rev(RColorBrewer::brewer.pal(3, "RdBu")),
                      # breaks = c(0, 1, max(res_mat_logp, na.rm = TRUE)),
                      color = rev(RColorBrewer::brewer.pal(11, "RdBu"))[c(1:5, 6, 8,11)],
                      breaks = c(seq(0, 1, by = 0.2), 1.01, max(res_mat_logp, na.rm = TRUE)),
                      cluster_rows = FALSE,
                      cluster_cols = FALSE))
  
  if(length(genes_names) > 50){
    pdf(file = paste0(out_dir, "pair_matches_summary_log.pdf"),
        onefile = FALSE,
        width = 18,
        height = 15)
  } else {
    pdf(file = paste0(out_dir, "pair_matches_summary_log.pdf"),
        onefile = FALSE,
        width = 12,
        height = 10)
  }
  print(p2)
  dev.off()
  
  # saveRDS(res_mat, file = paste0(out_dir, "pair_matches_mat.rds"))
  # saveRDS(res_df, file = paste0(out_dir, "pair_matches_df.rds"))
  
  colnames(res_df) <- c("gene_1", "gene_2", "p_value")
  res_df$p_value <- as.numeric(res_df$p_value)
  res_df <- res_df %>%
    as_tibble() %>%
    arrange(p_value)
  res_df_signif <- res_df %>%
    dplyr::filter(!is.na(p_value)) %>%
    dplyr::filter(p_value <= 0.05)
  res_df$p_value <- round(res_df$p_value, 4)
  res_df_signif$p_value <- round(res_df_signif$p_value, 4)
  write.csv(res_df, file = paste0(out_dir, "pair_matches_df.csv"))
  write.csv(res_df_signif, file = paste0(out_dir, "pair_matches_signif_df.csv"))
  write.csv(round(res_mat, 4), file = paste0(out_dir, "pair_matches_mat.csv"))
  
  
  
  
  # Save survival plots
  if(genes_type == "cis_top_hits"){
    res_df_nonSignif <- res_df %>%
      dplyr::filter(!is.na(p_value)) %>%
      dplyr::filter(p_value > 0.05)
    nonSignif_pairs <- paste0(res_df_nonSignif$gene_1, "&", res_df_nonSignif$gene_2)
    pdf(file = paste0(out_dir, "survival_nonSignif.pdf"),
        onefile = TRUE,
        width = 10,
        height = 12)
    print(nonSignif_plot_list[nonSignif_pairs])
    dev.off()
  }
  signif_pairs <- paste0(res_df_signif$gene_1, "&", res_df_signif$gene_2)
  pdf(file = paste0(out_dir, "survival_signif.pdf"),
      onefile = TRUE,
      width = 10,
      height = 12)
  print(signif_plot_list[signif_pairs])
  dev.off()
  
}

# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
# -----------------------------------------------------------------------------------
# ----------------------- SAVE TOP 15 GENE CIS INTEGRATED FUNCS ---------------------
# -----------------------------------------------------------------------------------
# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&

saveIntegratedGenesSurvival <- function(SU2CEC_metadata,
                                        SU2CEC_vpmat,
                                        genes_names,
                                        genes_vpscores,
                                        out_dir,
                                        out_name = "top_CIS_hits",
                                        thresh_method = "quantile"){
  # Add the integrated top CIS genes scores together
  n_samples <- length(SU2CEC_metadata$Sample_Identifier)
  genes_integrated_score <- rep(NA, n_samples)
  names(genes_integrated_score) <- SU2CEC_metadata$Sample_Identifier
  for(i in 1:n_samples){
    sample_i_name <- SU2CEC_metadata$Sample_Identifier[i]
    genes_integrated_score[i] <- sum(SU2CEC_vpmat[genes_names, sample_i_name]*sign(genes_vpscores))
  }
  SU2CEC_metadata <- cbind(SU2CEC_metadata, genes_integrated_score)
  colnames(SU2CEC_metadata)[ncol(SU2CEC_metadata)] <- out_name
  
  gene_i <- out_name
  
  if(thresh_method == "quantile"){
    SU2CEC_metadata_2 <- SU2CEC_metadata %>%
      dplyr::filter(get(gene_i) > quantile(SU2CEC_metadata[[gene_i]], 0.75) | get(gene_i) < quantile(SU2CEC_metadata[[gene_i]], 0.25))
  } else {
    SU2CEC_metadata_2 <- SU2CEC_metadata %>%
      dplyr::filter(get(gene_i) >= 1.96 | get(gene_i) <= -1.96)
  }
  
  SU2CEC_metadata_2[[gene_i]] <- as.integer(SU2CEC_metadata_2[[gene_i]] > 0)
  
  gene_i_coxph <- survival::coxph(Surv(time, status) ~ get(gene_i), data = SU2CEC_metadata_2)
  gene_i_coxph_plot <- gtsummary::tbl_regression(gene_i_coxph, exp = TRUE) 
  
  p.value <- gene_i_coxph_plot$table_body$p.value
  my_filename <- paste0(out_dir, gene_i, "_stats.txt")
  write(capture.output(print(gene_i_coxph)), my_filename)
  
  survival_plot <- survfit2(Surv(time, status) ~ get(gene_i), data = SU2CEC_metadata_2) %>% 
    ggsurvfit() +
    labs(
      x = "Days",
      y = "Overall survival probability"
    ) + 
    add_confidence_interval() +
    scale_color_manual(values = c("blue", "red"), labels=c(paste0(gene_i, "-"), paste0(gene_i, "+"))) +
    scale_fill_manual(values = c("blue", "red"), labels=c(paste0(gene_i, "-"), paste0(gene_i, "+"))) +
    add_risktable() +
    add_pvalue(location = "annotation", size = 5)
  
  my_filename <- paste0(out_dir, "/", gene_i, ".pdf")
  ggsave(filename = my_filename,
         plot = survival_plot,
         device = pdf(),
         onefile = FALSE)
  dev.off()
}

# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
# -----------------------------------------------------------------------------------
# --------------------------- NEPC VS ADENO SURVIVAL FUNC ---------------------------
# -----------------------------------------------------------------------------------
# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
saveNEPCvsAdenoSurvivalPlot <- function(SU2CEC_metadata, out_dir){
  survival_plot <- survfit2(Surv(time, status) ~ Neuroendocrine_Features, data = SU2CEC_metadata) %>% 
    ggsurvfit() +
    labs(
      x = "Days",
      y = "Overall survival probability"
    ) + 
    add_confidence_interval() +
    scale_color_manual(values = c("blue", "red"), labels=c("Adeno", "NEPC")) +
    scale_fill_manual(values = c("blue", "red"), labels=c("Adeno", "NEPC")) +
    add_risktable() +
    add_pvalue(location = "annotation", size = 5)
  my_filename <- paste0(out_dir, "NEPC_vs_Adeno_survival.pdf")
  ggsave(filename = my_filename,
         plot = survival_plot,
         device = pdf(),
         onefile = FALSE)
  dev.off()
}

# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
# -----------------------------------------------------------------------------------
# ------------------------ SU2C SURVIVAL ANALYSIS PIPELINE  -------------------------
# -----------------------------------------------------------------------------------
# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
performSurvivalAnalysisOnSU2CDataset <- function(SU2CEC_metadata,
                                                 SU2CEC_vpmat,
                                                 dataset_name = "SU2CEastCoast",
                                                 data_type = "pAct"){
  # Create folders
  reports_dir <- "experiments/survival-analysis/reports/"
  dir.create(reports_dir, showWarnings = FALSE)
  reports_dir <- paste0(reports_dir, dataset_name, "/")
  dir.create(reports_dir, showWarnings = FALSE)
  reports_dir <- paste0(reports_dir, data_type, "/")
  dir.create(reports_dir, showWarnings = FALSE)
  
  
  print("--------------------------------------------------------------------------")
  print("---------------------- Individual Survival Analysis ----------------------")
  print("--------------------------------------------------------------------------")
  dir.create(paste0(reports_dir, "individual_survival/"), showWarnings = FALSE)
  # Top CIS Genes
  dir.create(paste0(reports_dir, "individual_survival/top_cis_genes/"), showWarnings = FALSE)
  dir.create(paste0(reports_dir, "individual_survival/top_cis_genes/quantile_thresh/"), showWarnings = FALSE)
  getSurvivalForGeneList(SU2CEC_metadata,
                         SU2CEC_vpmat,
                         gene_list = cis_top_hits,
                         vpscores = cis_top_hits_vpscores,
                         out_dir = paste0(reports_dir, "individual_survival/top_cis_genes/quantile_thresh/"),
                         out_name = "top_CIS_genes",
                         thresh_method = "quantile")
  if(data_type == "pAct"){
    dir.create(paste0(reports_dir, "individual_survival/top_cis_genes/NES_thresh/"), showWarnings = FALSE)
    getSurvivalForGeneList(SU2CEC_metadata,
                           SU2CEC_vpmat,
                           gene_list = cis_top_hits,
                           vpscores = cis_top_hits_vpscores,
                           out_dir = paste0(reports_dir, "individual_survival/top_cis_genes/NES_thresh/"),
                           out_name = "top_CIS_genes",
                           thresh_method = "NES")
  }

  # Top 50 MRs
  dir.create(paste0(reports_dir, "individual_survival/top_50_mrs/"), showWarnings = FALSE)
  dir.create(paste0(reports_dir, "individual_survival/top_50_mrs/quantile_thresh/"), showWarnings = FALSE)
  getSurvivalForGeneList(SU2CEC_metadata,
                         SU2CEC_vpmat,
                         gene_list = top_candidate_mrs,
                         vpscores = top_candidate_mrs_vpscores,
                         out_dir = paste0(reports_dir, "individual_survival/top_50_mrs/quantile_thresh/"),
                         out_name = "top_50_mrs",
                         thresh_method = "quantile")
  if(data_type == "pAct"){
    dir.create(paste0(reports_dir, "individual_survival/top_50_mrs/NES_thresh/"), showWarnings = FALSE)
    getSurvivalForGeneList(SU2CEC_metadata,
                           SU2CEC_vpmat,
                           gene_list = top_candidate_mrs,
                           vpscores = top_candidate_mrs_vpscores,
                           out_dir = paste0(reports_dir, "individual_survival/top_50_mrs/NES_thresh/"),
                           out_name = "top_50_mrs",
                           thresh_method = "NES")
  }
  # Bottom 50 MRs
  dir.create(paste0(reports_dir, "individual_survival/bottom_50_mrs/"), showWarnings = FALSE)
  dir.create(paste0(reports_dir, "individual_survival/bottom_50_mrs/quantile_thresh/"), showWarnings = FALSE)
  getSurvivalForGeneList(SU2CEC_metadata,
                         SU2CEC_vpmat,
                         gene_list = bottom_candidate_mrs,
                         vpscores = bottom_candidate_mrs_vpscores,
                         out_dir = paste0(reports_dir, "individual_survival/bottom_50_mrs/quantile_thresh/"),
                         out_name = "bottom_50_mrs",
                         thresh_method = "quantile")
  if(data_type == "pAct"){
    dir.create(paste0(reports_dir, "individual_survival/bottom_50_mrs/NES_thresh/"), showWarnings = FALSE)
    getSurvivalForGeneList(SU2CEC_metadata,
                           SU2CEC_vpmat,
                           gene_list = bottom_candidate_mrs,
                           vpscores = bottom_candidate_mrs_vpscores,
                           out_dir = paste0(reports_dir, "individual_survival/bottom_50_mrs/NES_thresh/"),
                           out_name = "bottom_50_mrs",
                           thresh_method = "NES")
  }
  # All CIS Genes
  dir.create(paste0(reports_dir, "individual_survival/all_cis_genes/"), showWarnings = FALSE)
  dir.create(paste0(reports_dir, "individual_survival/all_cis_genes/quantile_thresh/"), showWarnings = FALSE)
  getSurvivalForGeneList(SU2CEC_metadata,
                         SU2CEC_vpmat,
                         gene_list = all_cis_genes,
                         vpscores = all_cis_genes_vpscores,
                         out_dir = paste0(reports_dir, "individual_survival/all_cis_genes/quantile_thresh/"),
                         out_name = "all_cis_genes",
                         thresh_method = "quantile")
  if(data_type == "pAct"){
    dir.create(paste0(reports_dir, "individual_survival/all_cis_genes/NES_thresh/"), showWarnings = FALSE)
    getSurvivalForGeneList(SU2CEC_metadata,
                           SU2CEC_vpmat,
                           gene_list = all_cis_genes,
                           vpscores = all_cis_genes_vpscores,
                           out_dir = paste0(reports_dir, "individual_survival/all_cis_genes/NES_thresh/"),
                           out_name = "all_cis_genes",
                           thresh_method = "NES")
  }

  print("--------------------------------------------------------------------------")
  print("----------------------- Pairwise Survival Analysis -----------------------")
  print("--------------------------------------------------------------------------")
  dir.create(paste0(reports_dir, "pairwise_survival"), showWarnings = FALSE)
  # Top CIS Genes
  dir.create(paste0(reports_dir, "pairwise_survival/pairwise_cis_genes/"), showWarnings = FALSE)
  dir.create(paste0(reports_dir, "pairwise_survival/pairwise_cis_genes/quantile_thresh/"), showWarnings = FALSE)
  performPairChecksAndSaveResults(SU2CEC_metadata = SU2CEC_metadata,
                                  SU2CEC_vpmat = SU2CEC_vpmat,
                                  genes_names = cis_top_hits,
                                  vpscores = cis_top_hits_vpscores,
                                  out_dir = paste0(reports_dir, "pairwise_survival/pairwise_cis_genes/quantile_thresh/"),
                                  genes_type = "cis_top_hits",
                                  thresh_method = "quantile")
  if(data_type == "pAct"){
    dir.create(paste0(reports_dir, "pairwise_survival/pairwise_cis_genes/NES_thresh/"), showWarnings = FALSE)
    performPairChecksAndSaveResults(SU2CEC_metadata = SU2CEC_metadata,
                                    SU2CEC_vpmat = SU2CEC_vpmat,
                                    genes_names = cis_top_hits,
                                    vpscores = cis_top_hits_vpscores,
                                    out_dir = paste0(reports_dir, "pairwise_survival/pairwise_cis_genes/NES_thresh/"),
                                    genes_type = "cis_top_hits",
                                    thresh_method = "NES")
  }
  # Top MRs
  dir.create(paste0(reports_dir, "pairwise_survival/pairwise_mrs/"), showWarnings = FALSE)
  dir.create(paste0(reports_dir, "pairwise_survival/pairwise_mrs/quantile_thresh/"), showWarnings = FALSE)
  performPairChecksAndSaveResults(SU2CEC_metadata = SU2CEC_metadata,
                                  SU2CEC_vpmat = SU2CEC_vpmat,
                                  genes_names = all_mrs_names,
                                  vpscores = all_mrs_vpscores,
                                  out_dir = paste0(reports_dir, "pairwise_survival/pairwise_mrs/quantile_thresh/"),
                                  save_plots = FALSE,
                                  genes_type = "top_mrs",
                                  thresh_method = "quantile")
  if(data_type == "pAct"){
    dir.create(paste0(reports_dir, "pairwise_survival/pairwise_mrs/NES_thresh/"), showWarnings = FALSE)
    performPairChecksAndSaveResults(SU2CEC_metadata = SU2CEC_metadata,
                                    SU2CEC_vpmat = SU2CEC_vpmat,
                                    genes_names = all_mrs_names,
                                    vpscores = all_mrs_vpscores,
                                    out_dir = paste0(reports_dir, "pairwise_survival/pairwise_mrs/NES_thresh/"),
                                    save_plots = FALSE,
                                    genes_type = "top_mrs",
                                    thresh_method = "NES")
  }

  print("--------------------------------------------------------------------------")
  print("-------------------- Top 15 Genes Integrated Analysis --------------------")
  print("--------------------------------------------------------------------------")
  dir.create(paste0(reports_dir, "integrated_survival/"), showWarnings = FALSE)
  dir.create(paste0(reports_dir, "integrated_survival/quantile_thresh/"), showWarnings = FALSE)
  saveIntegratedGenesSurvival(SU2CEC_metadata,
                              SU2CEC_vpmat,
                              genes_names = cis_top_hits,
                              genes_vpscores = cis_top_hits_vpscores,
                              out_dir = paste0(reports_dir, "integrated_survival/quantile_thresh/"),
                              out_name = "top_15_CIS_hits",
                              thresh_method = "quantile")
  if(data_type == "pAct"){
    dir.create(paste0(reports_dir, "integrated_survival/NES_thresh/"), showWarnings = FALSE)
    saveIntegratedGenesSurvival(SU2CEC_metadata,
                             SU2CEC_vpmat,
                             genes_names = cis_top_hits,
                             genes_vpscores = cis_top_hits_vpscores,
                             out_dir = paste0(reports_dir, "integrated_survival/NES_thresh/"),
                             out_name = "top_15_CIS_hits",
                             thresh_method = "NES")
  }

  print("--------------------------------------------------------------------------")
  print("----------- Top 6-7 POSITIVE TOP CIS Genes Integrated Analysis -----------")
  print("--------------------------------------------------------------------------")
  dir.create(paste0(reports_dir, "integrated_pos_survival/"), showWarnings = FALSE)
  dir.create(paste0(reports_dir, "integrated_pos_survival/quantile_thresh/"), showWarnings = FALSE)
  cis_top_pos_hits <- names(sort(cis_top_hits_vpscores, decreasing = TRUE)[1:6])
  saveIntegratedGenesSurvival(SU2CEC_metadata,
                              SU2CEC_vpmat,
                              genes_names = cis_top_pos_hits,
                              genes_vpscores = cis_top_hits_vpscores[cis_top_pos_hits],
                              out_dir = paste0(reports_dir, "integrated_pos_survival/quantile_thresh/"),
                              out_name = "top_pos_CIS_hits",
                              thresh_method = "quantile")
  if(data_type == "pAct"){
    dir.create(paste0(reports_dir, "integrated_pos_survival/NES_thresh/"), showWarnings = FALSE)
    saveIntegratedGenesSurvival(SU2CEC_metadata,
                                SU2CEC_vpmat,
                                genes_names = cis_top_pos_hits,
                                genes_vpscores = cis_top_hits_vpscores[cis_top_pos_hits],
                                out_dir = paste0(reports_dir, "integrated_pos_survival/NES_thresh/"),
                                out_name = "top_pos_CIS_hits",
                                thresh_method = "NES")
  }



  print("--------------------------------------------------------------------------")
  print("-------------------------------- NEPC Plot -------------------------------")
  print("--------------------------------------------------------------------------")
  dir.create(paste0(reports_dir, "phenotype_survival/"), showWarnings = FALSE)
  saveNEPCvsAdenoSurvivalPlot(SU2CEC_metadata, out_dir = paste0(reports_dir, "phenotype_survival/"))
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

# Load VIPER Data
# Beltran_vpmat <- readRDS("experiments/oncomatch-analysis/processed_data/viper/Beltran_viper_zscore_SU2CNets.rds")
# SB_vpmat <- readRDS("experiments/oncomatch-analysis/processed_data/viper/SleepingBeauty_viper_zscore_SU2CNets_SBNet.rds")
SU2CEC_vpmat <- readRDS("experiments/oncomatch-analysis/processed_data/viper/SU2CEastCoast_viper_zscore_SU2CNets.rds")
SU2CWC_vpmat <- readRDS("experiments/oncomatch-analysis/processed_data/viper/SU2CWestCoast_viper_zscore_SU2CNets.rds")
SU2CEC_gExpr <- GESTransform(readRDS("experiments/oncomatch-analysis/processed_data/gexpr/SU2CEastCoast_gExpr_polyA_tpm.rds"))
SU2CWC_gExpr <- GESTransform(readRDS("data/SU2C_WestCoast/SU2CWestCoast_merged_tpm_geneNames.rds"))
colnames(SU2CWC_vpmat) = colnames(SU2CWC_vpmat) %>%
  stringr::str_replace_all("_2018", "") %>%
  stringr::str_replace_all("\\.", "-")
colnames(SU2CWC_gExpr) = colnames(SU2CWC_gExpr) %>%
  stringr::str_replace_all("_2018", "") %>%
  stringr::str_replace_all("\\.", "-")

# Get SU2C EC and WC Metadata
SU2CEC_metadata <- loadSu2CECMetadata()
SU2CWC_metadata <- loadSu2CWCMetadata()

colnames(SU2CEC_gExpr)[colnames(SU2CEC_vpmat) %in% SU2CEC_metadata$Sample_Identifier] <- SU2CEC_metadata$Sample_Identifier

# Load CIS Top Hits And VPScores
my_top_hits_table <- readxl::read_excel("experiments/integrative-analysis/reports/nepc-signatures-table-all-3-conditions.xlsx")
my_cis_top_hits_list <- getCisTopHitsAndVPScores(my_top_hits_table)
cis_top_hits <- my_cis_top_hits_list[["cis_top_hits"]]
cis_top_hits_vpscores <- my_cis_top_hits_list[["cis_top_hits_vpscores"]]

# Load Top & Bottom 50 MRs And Signs (as VPScores)
my_table <- getMyTable("experiments/integrative-analysis/processed_data/tibble-of-sb-cis-nepc-viper-cindy-integrative-analysis.rds")
# Top MRs
top_candidate_mrs <- getTFsWithHighestNEPCScoresDF(my_table, n_top = 50)$human_gene_name
top_candidate_mrs_vpscores <- rep(1, 50)
names(top_candidate_mrs_vpscores) <- top_candidate_mrs
# Bottom MRs
bottom_candidate_mrs_vpscores <- rep(-1, 50)
bottom_candidate_mrs <- getTFsWithLowestNEPCScoresDF(my_table, n_lowest = 50)$human_gene_name
names(bottom_candidate_mrs_vpscores) <- bottom_candidate_mrs
# Combined MRs
all_mrs_names <- c(top_candidate_mrs, bottom_candidate_mrs)
all_mrs_vpscores <- c(rep(1, 50), rep(-1, 50))
names(all_mrs_vpscores) <- all_mrs_names
# All CIS Genes
all_cis_genes <- my_table %>%
  filter(is_cis_gene == "Yes") %>%
  pull(human_gene_name) %>%
  na.omit() %>%
  as.character()
all_cis_genes_vpscores <- my_table %>%
  filter(is_cis_gene == "Yes") %>%
  pull(NEPC_viper_score) 
names(all_cis_genes_vpscores) <- all_cis_genes
all_cis_genes_vpscores[is.na(all_cis_genes_vpscores)] <- 1

# Perform Survival Analysis
performSurvivalAnalysisOnSU2CDataset(SU2CEC_metadata, SU2CEC_vpmat, dataset_name = "SU2CEastCoast", data_type = "pAct")
performSurvivalAnalysisOnSU2CDataset(SU2CWC_metadata, SU2CWC_vpmat, dataset_name = "SU2CWestCoast", data_type = "pAct")

performSurvivalAnalysisOnSU2CDataset(SU2CEC_metadata, SU2CEC_gExpr, dataset_name = "SU2CEastCoast", data_type = "gExpr")
performSurvivalAnalysisOnSU2CDataset(SU2CWC_metadata, SU2CWC_gExpr, dataset_name = "SU2CWestCoast", data_type = "gExpr")

SU2C_EC_gExpr_quantThresh_pvals <- read.csv("experiments/survival-analysis/reports/SU2CEastCoast/gExpr/individual_survival/all_cis_genes/quantile_thresh/all_cis_genes_pvals_df.csv", row.names = 1) %>%
  rbind(
    read.csv("experiments/survival-analysis/reports/SU2CEastCoast/gExpr/individual_survival/top_50_mrs/quantile_thresh/top_50_mrs_pvals_df.csv", row.names = 1)
  ) %>%
  rbind(
    read.csv("experiments/survival-analysis/reports/SU2CEastCoast/gExpr/individual_survival/bottom_50_mrs/quantile_thresh/bottom_50_mrs_pvals_df.csv", row.names = 1)
  ) %>%
  as_tibble() %>%
  dplyr::rename("SU2C_EC_gExpr_quantThresh_survival_CIS_pvals" = "p_value")
SU2C_EC_pAct_quantThresh_pvals <- read.csv("experiments/survival-analysis/reports/SU2CEastCoast/pAct/individual_survival/all_cis_genes/quantile_thresh/all_cis_genes_pvals_df.csv", row.names = 1) %>%
  rbind(
    read.csv("experiments/survival-analysis/reports/SU2CEastCoast/pAct/individual_survival/top_50_mrs/quantile_thresh/top_50_mrs_pvals_df.csv", row.names = 1)
  ) %>%
  rbind(
    read.csv("experiments/survival-analysis/reports/SU2CEastCoast/pAct/individual_survival/bottom_50_mrs/quantile_thresh/bottom_50_mrs_pvals_df.csv", row.names = 1)
  ) %>%
  as_tibble() %>%
  dplyr::rename("SU2C_EC_pAct_quantThresh_survival_CIS_pvals" = "p_value")
SU2C_EC_pAct_NESThresh_pvals <- read.csv("experiments/survival-analysis/reports/SU2CEastCoast/pAct/individual_survival/all_cis_genes/NES_thresh/all_cis_genes_pvals_df.csv", row.names = 1) %>%
  rbind(
    read.csv("experiments/survival-analysis/reports/SU2CEastCoast/pAct/individual_survival/top_50_mrs/NES_thresh/top_50_mrs_pvals_df.csv", row.names = 1)
  ) %>%
  rbind(
    read.csv("experiments/survival-analysis/reports/SU2CEastCoast/pAct/individual_survival/bottom_50_mrs/NES_thresh/bottom_50_mrs_pvals_df.csv", row.names = 1)
  ) %>%
  as_tibble() %>%
  dplyr::rename("SU2C_EC_pAct_NESThresh_survival_CIS_pvals" = "p_value")
SU2C_WC_gExpr_quantThresh_pvals <- read.csv("experiments/survival-analysis/reports/SU2CWestCoast/gExpr/individual_survival/all_cis_genes/quantile_thresh/all_cis_genes_pvals_df.csv", row.names = 1) %>%
  rbind(
    read.csv("experiments/survival-analysis/reports/SU2CWestCoast/gExpr/individual_survival/top_50_mrs/quantile_thresh/top_50_mrs_pvals_df.csv", row.names = 1)
  ) %>%
  rbind(
    read.csv("experiments/survival-analysis/reports/SU2CWestCoast/gExpr/individual_survival/bottom_50_mrs/quantile_thresh/bottom_50_mrs_pvals_df.csv", row.names = 1)
  ) %>%
  as_tibble() %>%
  dplyr::rename("SU2C_WC_gExpr_quantThresh_survival_CIS_pvals" = "p_value")
SU2C_WC_pAct_quantThresh_pvals <- read.csv("experiments/survival-analysis/reports/SU2CWestCoast/pAct/individual_survival/all_cis_genes/quantile_thresh/all_cis_genes_pvals_df.csv", row.names = 1) %>%
  rbind(
    read.csv("experiments/survival-analysis/reports/SU2CWestCoast/pAct/individual_survival/top_50_mrs/quantile_thresh/top_50_mrs_pvals_df.csv", row.names = 1)
  ) %>%
  rbind(
    read.csv("experiments/survival-analysis/reports/SU2CWestCoast/pAct/individual_survival/bottom_50_mrs/quantile_thresh/bottom_50_mrs_pvals_df.csv", row.names = 1)
  ) %>%
  as_tibble() %>%
  dplyr::rename("SU2C_WC_pAct_quantThresh_survival_CIS_pvals" = "p_value")
SU2C_WC_pAct_NESThresh_pvals <- read.csv("experiments/survival-analysis/reports/SU2CWestCoast/pAct/individual_survival/all_cis_genes/NES_thresh/all_cis_genes_pvals_df.csv", row.names = 1) %>%
  rbind(
    read.csv("experiments/survival-analysis/reports/SU2CWestCoast/pAct/individual_survival/top_50_mrs/NES_thresh/top_50_mrs_pvals_df.csv", row.names = 1)
  ) %>%
  rbind(
    read.csv("experiments/survival-analysis/reports/SU2CWestCoast/pAct/individual_survival/bottom_50_mrs/NES_thresh/bottom_50_mrs_pvals_df.csv", row.names = 1)
  ) %>%
  as_tibble() %>%
  dplyr::rename("SU2C_WC_pAct_NESThresh_survival_CIS_pvals" = "p_value")

my_table_updated <- my_table %>%
  dplyr::left_join(SU2C_EC_gExpr_quantThresh_pvals, by = c("human_gene_name" = "gene")) %>%
  dplyr::left_join(SU2C_EC_pAct_quantThresh_pvals, by = c("human_gene_name" = "gene")) %>%
  dplyr::left_join(SU2C_EC_pAct_NESThresh_pvals, by = c("human_gene_name" = "gene")) %>%
  dplyr::left_join(SU2C_WC_gExpr_quantThresh_pvals, by = c("human_gene_name" = "gene")) %>%
  dplyr::left_join(SU2C_WC_pAct_quantThresh_pvals, by = c("human_gene_name" = "gene")) %>%
  dplyr::left_join(SU2C_WC_pAct_NESThresh_pvals, by = c("human_gene_name" = "gene"))

writexl::write_xlsx(my_table_updated, "experiments/survival-analysis/reports/tibble-of-sb-cis-nepc-viper-cindy-survival-integrative-analysis.xlsx")

print("Done.")
