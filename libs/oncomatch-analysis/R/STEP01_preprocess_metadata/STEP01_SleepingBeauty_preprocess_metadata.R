suppressPackageStartupMessages(library(tidyverse))

# gemm_metadata <- readRDS("data/oncomatch-analysis/SleepingBeauty/sb-rna-seq-data-analysis.rds")$metadata %>%
gemm_metadata <- readRDS("experiments/all-samples-sb-rna-seq-analysis/processed_data/sb-rna-seq-data-analysis.rds")$metadata %>%
  as_tibble() %>%
  dplyr::filter(genotype == "NPp53") %>%
  # dplyr::select(-c("beltran_enrichment", "novel_beltran_enrichment")) %>%
  dplyr::select(-c("beltran_enrichment")) %>%
  dplyr::rename("Sample_Identifier" = "sample_id",
                "Phenotype" = "phenotype",
                "SB_Status" = "sb_status") 

saveRDS(gemm_metadata, "experiments/oncomatch-analysis/processed_data/metadata/SleepingBeauty_metadata_forOMPlot.rds")
