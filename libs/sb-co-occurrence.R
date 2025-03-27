reports.dir <- "experiments/table"

library(data.table)
library(stringr)

# Modified somatic interactions function
source("libs/tools/somaticInteractions.R")

# Function to count the number of entries in a row
countentries <- function(x) {sum(!grepl("[0]", x))}

# Function to count the number of altered samples per gene
countgenes <- function(df){
  tempvec <- c()
  for(i in 1:nrow(df)){
    tempvec[i] <- as.numeric(countentries(df[i,]))
  }
  genetotals <- cbind(rownames(df), as.numeric(tempvec))
  colnames(genetotals) <- c("Hugo_Symbol", "AlteredSamples")
  return(data.table(genetotals))
}

# Read in the metadata and full TAPDANCE output file
metadata <- read.csv("data/tapdance/mouse_metadata_20220912_FA.csv") 
mutation_raw <- read.csv("experiments/table/TAPDANCE_001.csv")

# Function to prepare the mutation matrix 
mutation_matrix_prep <- function(genes){
  mutation_red <- mutation_raw %>% dplyr::filter(gene_name %in% genes)
  
  mutation_df <- as.data.frame(matrix(0, nrow = length(genes), ncol = nrow(metadata)))
  rownames(mutation_df) <- genes
  colnames(mutation_df) <- metadata$DNA.Sample.ID[metadata$DNA.Sample.ID != ""]
  
  for(i in 1:nrow(mutation_red)){
    gene <- mutation_red$gene_name[i]
    samples <- trimws(unlist(strsplit(mutation_red$library_name[i], "::")), which = "both")
    samples <- samples[samples != ""]
    for(j in 1:ncol(mutation_df)){
      if(colnames(mutation_df)[j] %in% samples){mutation_df[gene,j] <- 1}
    }
  }
  mutation_df[is.na(mutation_df)] <- 0
  
  total_counts_sb <<- countgenes(mutation_df)
  total_counts_sb$AlteredSamples <<- as.numeric(total_counts_sb$AlteredSamples)
  genes_sb <<- total_counts_sb %>% dplyr::arrange(desc(AlteredSamples)) %>% dplyr::pull(Hugo_Symbol)
  mutation_matrix_sb <<- mutation_df[rownames(mutation_df) %in% genes_sb,]
  all_cases_sb <<- metadata$DNA.Sample.ID
  
}

# Get the genes that are in the CIS
CIS_genes <- read.csv("experiments/table/TAPDANCE_001.csv") %>% dplyr::distinct(pos, .keep_all = TRUE) %>% dplyr::pull(gene_name)

# Prepare the mutation matrix
mutation_matrix_prep(CIS_genes)

# Run the somatic interactions function
somatic_sb <- somaticInteractions(oncomatrix_passed = mutation_matrix_sb, 
                                  all.tsbs_passed = all_cases_sb,
                                  gene_sum_passed = total_counts_sb,
                                  genes = genes_sb,
                                  plot = FALSE,
                                  fontSize = .5)

write.csv(somatic_sb, file.path(reports.dir,"SB_co-occurrence_CIS_genes.csv"))

