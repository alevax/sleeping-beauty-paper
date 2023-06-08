suppressPackageStartupMessages(library(tidyverse))

beltran_clinical <- readRDS("data/Beltran-2016/beltran-clinical.rds")
beltran_clinical <- beltran_clinical %>%
  dplyr::select("sample_id",
                "Tumor Disease Anatomic Site",
                "Disease code",
                "Pathology Classification") %>% 
  rename("Sample_Identifier" = "sample_id",
         "Neuroendocrine_Features" = "Disease code",
         "Pathology_Classification" = "Pathology Classification",
         "Tissue_Site" = "Tumor Disease Anatomic Site") %>%
  mutate(Neuroendocrine_Features = case_when(
    Neuroendocrine_Features %in% c("CRPC-NE") ~ "Yes",
    Neuroendocrine_Features %in% c("CRPC-Adeno") ~ "No"
  )) %>%
  mutate(Pathology_Classification = case_when(
    Pathology_Classification == "A" ~ "Adenocarcinoma",
    Pathology_Classification == "B" ~ "Adenocarcinoma with NE features",
    Pathology_Classification == "C" ~ "Pure small cell",
    Pathology_Classification == "D" ~ "Large cell NE carcinoma",
    Pathology_Classification == "E" ~ "Mixed small cell"
  ))
# Beltran, H., Prandi, D., Mosquera, J. M., Benelli, M., Puca, L., Cyrta, J., ... & Demichelis, F. (2016). Divergent clonal evolution of castration-resistant neuroendocrine prostate cancer. Nature medicine, 22(3), 298-305.
  # Category A represents usual prostate adenocarcinoma without neuroendocrine differentiation
  # Category B represents usual prostate adenocarcinoma with neuroendocrine differentiation > 20%
  # Category C represents small cell carcinoma
  # Category D represents large cell neuroendocrine carcinoma
  # Category E represents mixed small cell carcinoma – adenocarcinoma.

beltran_clinical$Tissue_Site <- beltran_clinical$Tissue_Site %>%
  stringr::str_replace_all("Lymph node", "LN") %>%
  stringr::str_replace_all("Soft tissue", "Other Soft tissue")
beltran_clinical$NEPC_Score <- rep(NA, nrow(beltran_clinical))

saveRDS(beltran_clinical, "experiments/oncomatch-analysis/processed_data/metadata/Beltran_metadata_forOMPlot.rds")
