reports.dir <- "experiments/table"

library(tidyverse)
library(scales)

# Function to split genes into separate rows
split_genes <- function(df){
  td_df <- df %>% mutate_if(is.factor, as.character)
  td_df_clean <- data.frame(matrix(ncol = length(colnames(td_df))))
  colnames(td_df_clean) <- colnames(td_df)
  for(i in 1:length(td_df$gene_name)) {
    genes <- unlist(strsplit(td_df$gene_name[i], ","))
    for (j in 1:length(genes)){
      if (!grepl("^\\s*$", genes[j])) {
        rowtoadd <- td_df[i,]
        rowtoadd["gene_name"] <- genes[j]
        td_df_clean <- rbind(td_df_clean, rowtoadd)
      }
    }
  }
  td_df_clean <- td_df_clean[-1,]
  td_df_clean[,"gene_name"] <- td_df_clean[,"gene_name"] %>% trimws()
  return(td_df_clean)
}

# Function to add metadata to the split table
metadata_TD <- function(df, split_string = "::"){
  df$NE <- NULL
  df$AD <- NULL
  df$NE_proportion <- NULL
  df$NE_p.value <- NULL
 
  for(i in 1:nrow(df)){
    df$CIS_chr[i] <- strsplit(df$pos[i], ":")[[1]][1]
    df$CIS_start[i] <- strsplit(strsplit(df$pos[i], ":")[[1]][2],"-")[[1]][1]
    df$CIS_end[i] <- strsplit(strsplit(df$pos[i], ":")[[1]][2],"-")[[1]][2]
    samples <- unlist(strsplit(df$library_name[i], split_string))
    df$NE[i] <- sum(samples %in% NE_ids) 
    df$AD[i] <- sum(samples %in% AD_ids)
    df$NE_proportion[i] <- binom.test(x = df$NE[i], n = df$X.library[i], p = p_NE)$estimate
    df$NE_p.value[i] <- binom.test(x = df$NE[i], n = df$X.library[i], p = p_NE)$p.value
    df$NE_prediction <- ifelse(df$NE_p.value < .05, ifelse(df$NE_proportion > p_NE,"Neuroendocrine","Adenocarcinoma"), NA)
    df$Drive_proportion[i] <- binom.test(x = df$number.of.inserts.drive.transcription.on.positive.strand[i], n = df$X.inserts[i], p = p_NE)$estimate
    df$Drive_p.value[i] <- binom.test(x = df$number.of.inserts.drive.transcription.on.positive.strand[i], n = df$X.inserts[i], p = p_NE)$p.value
    df$Drive_prediction <- ifelse(df$Drive_p.value < .05, ifelse(df$Drive_proportion > .5,"Driving","Disrupting"), NA)
  }
  df %>% mutate(NE_proportion=try(percent(NE_proportion),silent = T),
                Drive_proportion=try(percent(Drive_proportion),silent = T)) %>% return()
}

# Function to add metadata to the unsplit table
add_metadata_to_unsplit <- function(df, split_string = "::") {
  df$NE <- 0
  df$AD <- 0
  df$NE_proportion <- NA
  df$NE_p.value <- NA
  df$NE_prediction <- NA
  df$Drive_proportion <- NA
  df$Drive_p.value <- NA
  df$Drive_prediction <- NA
  df$CIS_chr <- NA
  df$CIS_start <- NA
  df$CIS_end <- NA
  
  for (i in 1:nrow(df)) {
    # Split CIS position
    pos_split <- strsplit(df$pos[i], ":")[[1]]
    df$CIS_chr[i] <- pos_split[1]
    range_split <- strsplit(pos_split[2], "-")[[1]]
    df$CIS_start[i] <- range_split[1]
    df$CIS_end[i] <- range_split[2]
    
    # Count NE and AD samples
    samples <- unlist(strsplit(df$library_name[i], split_string))
    df$NE[i] <- sum(samples %in% NE_ids)
    df$AD[i] <- sum(samples %in% AD_ids)
    
    # Calculate NE proportion and p-value
    ne_test <- binom.test(x = df$NE[i], n = df$X.library[i], p = p_NE)
    df$NE_proportion[i] <- ne_test$estimate
    df$NE_p.value[i] <- ne_test$p.value
    
    # Determine NE prediction
    df$NE_prediction[i] <- ifelse(df$NE_p.value[i] < .05, 
                                  ifelse(df$NE_proportion[i] > p_NE, "Neuroendocrine", "Adenocarcinoma"), 
                                  NA)
    
    # Calculate Drive proportion and p-value
    drive_test <- binom.test(x = df$number.of.inserts.drive.transcription.on.positive.strand[i], 
                             n = df$X.inserts[i], 
                             p = 0.5)
    df$Drive_proportion[i] <- drive_test$estimate
    df$Drive_p.value[i] <- drive_test$p.value
    
    # Determine Drive prediction
    df$Drive_prediction[i] <- ifelse(df$Drive_p.value[i] < .05, 
                                     ifelse(df$Drive_proportion[i] > 0.5, "Driving", "Disrupting"), 
                                     NA)
  }
  
  # Format proportions as percentages
  df$NE_proportion <- scales::percent(df$NE_proportion)
  df$Drive_proportion <- scales::percent(df$Drive_proportion)
  
  return(df)
}

# Read in file
TD.001 <- read.csv("TAPDANCE_NEPC/results/all/cis_all-nr-prostate-0.001.csv", header = T)

# Split into different genes per line
TD.001_split <- split_genes(TD.001) %>% unique()

# Look in metadata to see mouse phenotype
metadata <- read.csv("data/tapdance/mouse_metadata_20220912_FA.csv", header = TRUE)
NE_ids <- metadata %>% dplyr::filter(NE.Phenotype == TRUE) %>% pull(DNA.Sample.ID)
AD_ids <- metadata %>% dplyr::filter(NE.Phenotype == FALSE) %>% pull(DNA.Sample.ID)

# Determine baseline proportion of NE samples
p_NE <- sum(metadata$NE.Phenotype, na.rm = TRUE)/(nrow(metadata)-sum(is.na(metadata$NE.Phenotype)))

TD.001_metadata <- TD.001_split %>% metadata_TD()
write.csv(TD.001_metadata, file.path(reports.dir,"TAPDANCE_001.csv"))

TD.001_metadata_unsplit <- add_metadata_to_unsplit(TD.001)
write.csv(TD.001_metadata_unsplit, file.path(reports.dir, "TAPDANCE_001_unsplit.csv"))