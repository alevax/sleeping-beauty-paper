# Libraries
suppressPackageStartupMessages(library(tidyverse))

# ------------------------------------------------------------------------------------------
# -------------------------------- GET THE CLINICAL METADATA -------------------------------
# ------------------------------------------------------------------------------------------
# Load in the oncomatch-analysis metadata
european_urology_2019 <- as_tibble(read.csv("data/SU2C_WestCoast/WCDT/clinical_metadata/2021_05_13_european_urology_2019.txt", sep = "\t"))
JCO_2018 <- as_tibble(read.csv("data/SU2C_WestCoast/WCDT/clinical_metadata/2021_11_21_JCO_2018.txt", sep = "\t"))
JCO_tissue_sites <- as_tibble(read.csv("data/SU2C_WestCoast/WCDT/clinical_metadata/2022_10_18_JCO_tissue_sites.txt", sep = "\t"))

# Load in the processed merged counts
merged_counts_colnames <- colnames(readRDS("data/SU2C_WestCoast/SU2CWestCoast_merged_tpm_geneNames.rds"))
merged_counts_colnames <- merged_counts_colnames %>%
  stringr::str_replace_all("_2018", "") %>%
  stringr::str_replace_all("\\.", "-")

# Update sample IDs to keep them consistent
JCO_2018$sample_id[JCO_2018$sample_id == "DTB-098-PRO"] <- "DTB-098-PRO2"
JCO_2018$sample_id[JCO_2018$sample_id == "DTB-092-BL"] <- "DTB-092-BL-2018"
JCO_2018$sample_id[JCO_2018$sample_id == "DTB-055-PRO2"] <- "DTB-055-PRO"
european_urology_2019$sample_id[european_urology_2019$sample_id == "DTB-092-BL"] <- "DTB-092-BL-2018"

# Filter sample IDs to those in the merged_counts
JCO_2018 <- JCO_2018 %>%
  filter(sample_id %in% merged_counts_colnames)
european_urology_2019 <- european_urology_2019 %>%
  filter(sample_id %in% merged_counts_colnames)

clinical_metadata <- JCO_2018 %>%
  left_join(JCO_tissue_sites, by = c("sample_id" = "sample_id")) %>%
  left_join(european_urology_2019, by = c("sample_id" = "sample_id"))
clinical_metadata$id_aggarwal <- NULL

# ------------------------------------------------------------------------------------------
# --------------------- RESTRUCTURE THE CLINICAL METADATA FOR PLOTTING ---------------------
# ------------------------------------------------------------------------------------------
# print("Add pathology call to the metadata")
cell_pathology_metadata <- as_tibble(read.csv("data/SU2C_WestCoast/cell_pathology_SU2C_West_Coast.csv"))
clinical_metadata <- dplyr::left_join(clinical_metadata, cell_pathology_metadata, by = c("sample_id" = "sample_id"))

clinical_metadata_sel <- clinical_metadata %>%
  dplyr::select(sample_id,
                cluster_ID,
                OS_duration,
                "OS_event_.dead.",
                prior_Abi,
                prior_Enza,
                enzalutamide.resistant,
                metastasis_biopsy_site,
                pathology_call
  ) %>%
  mutate(Neuroendocrine_Features = case_when(
    pathology_call %in% c("Pure small cell", "Mixed small cell") ~ TRUE,#"Yes",
    pathology_call %in% c("Not small cell") ~ FALSE,#"No",
    pathology_call %in% c("Unclassified") ~ NA
  )) %>%
  mutate(Neuroendocrine_Features = as.character(Neuroendocrine_Features)) %>%
  mutate(Neuroendocrine_Features = stringr::str_replace_all(Neuroendocrine_Features,
                                                            pattern = "TRUE",
                                                            replacement = "Yes")) %>%
  mutate(Neuroendocrine_Features = stringr::str_replace_all(Neuroendocrine_Features,
                                                            pattern = "FALSE",
                                                            replacement = "No")) %>%
  rename("Overall_Survival_Days" = "OS_duration",
         "Overall_Survival_Status" = "OS_event_.dead.",
         "Pathology_Classification" = "pathology_call",
         "Tissue_Site" = "metastasis_biopsy_site",
         "Sample_Identifier" = "sample_id")

clinical_metadata_sel <- clinical_metadata_sel %>%
  mutate(Drug_Exposure_Status = case_when(
    prior_Abi == "Resistant" | prior_Enza == "Resistant" ~ "Resistant",
    prior_Abi == "Naïve" & prior_Enza == "Naïve" ~ "Naïve"
  ))

clinical_metadata_sel$Tissue_Site <- clinical_metadata_sel$Tissue_Site %>%
  stringr::str_to_title() %>%
  stringr::str_replace_all("Other_soft_tissue", "Other Soft tissue") %>%
  stringr::str_replace_all("Ln", "LN")

clinical_metadata_sel$Pathology_Classification <- clinical_metadata_sel$Pathology_Classification %>%
  stringr::str_replace_all("Not small cell", "Adenocarcinoma")

su2cWC_metadata_toPlot <- clinical_metadata_sel %>%
  dplyr::select("Sample_Identifier",
                "Tissue_Site",
                "Pathology_Classification",
                "Neuroendocrine_Features"
  )

su2cWC_metadata_toPlot$NEPC_Score <- rep(NA, nrow(su2cWC_metadata_toPlot))

saveRDS(su2cWC_metadata_toPlot, "experiments/oncomatch-analysis/processed_data/metadata/SU2C_WestCoast_metadata_forOMPlot.rds")
