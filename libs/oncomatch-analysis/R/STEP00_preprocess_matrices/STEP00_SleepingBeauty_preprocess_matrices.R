# Libraries
suppressPackageStartupMessages(library(tidyverse))

# Load files
gemm_samples <- readRDS("data/sb-counts/gemms-all-samples-tpm.rds")
mouse2human <- as_tibble(read.csv("data/translate/mouse2human.csv", row.names = 1))
indices_to_delete <- c()

pb = txtProgressBar(min = 0, max = nrow(gemm_samples), initial = 0, style = 3)
keep_row <- rep(TRUE, nrow(gemm_samples))
for(i in 1:nrow(gemm_samples)){
  gene_i_mouse_symbol <- rownames(gemm_samples)[i]
  
  gene_i_human_symbol <- mouse2human %>%
    filter(mouse_symbol == gene_i_mouse_symbol) %>%
    slice(1) %>%
    pull(human_symbol)
  if(length(gene_i_human_symbol) == 0){
    keep_row[i] <- FALSE
  } else {
    rownames(gemm_samples)[i] <- gene_i_human_symbol
  }
  setTxtProgressBar(pb, i)
}
close(pb)
gemm_samples <- gemm_samples[keep_row, ]
# dim(gemm_samples)
# [1] 16833   201
saveRDS(gemm_samples, "experiments/oncomatch-analysis/processed_data/gexpr/SleepingBeauty_gExpr_tpm_humanNames.rds")
