suppressPackageStartupMessages(library(tidyverse))

# print("Loading data files...")
data_clinical_patient <- as_tibble(read.csv("data/SU2C_EastCoast/prad_su2c_2019/data_clinical_patient.txt", sep = "\t"))
data_clinical_sample <- as_tibble(read.csv("data/SU2C_EastCoast/prad_su2c_2019/data_clinical_sample.txt", sep = "\t"))
gExpr <- readRDS("experiments/oncomatch-analysis/processed_data/gexpr/SU2CEastCoast_gExpr_polyA_tpm.rds")

# print("Reformatting the clinical patient data...")
data_clinical_patient <- data_clinical_patient %>%
  slice(-c(1:3)) %>%
  dplyr::rename("Overall_Survival_Months" = "Overall.Survival..Months.",
                "Chemo_Regimen_Category" = "Chemo.Regimen.Category",
                "Patient_Identifier" = "X.Patient.Identifier",
                "Overall_Survival_Status" = "Overall.Survival.Status") %>%
  mutate(Overall_Survival_Months = suppressWarnings(as.numeric(Overall_Survival_Months))) %>%
  mutate(Overall_Survival_Days = Overall_Survival_Months*30.4167) %>%
  dplyr::select(
    Patient_Identifier,
    Chemo_Regimen_Category,
    Overall_Survival_Status,
    Overall_Survival_Days
  ) %>%
  slice(-1)

# print("Reformatting the clinical sample data...")
data_clinical_sample <- data_clinical_sample %>%
  slice(-c(1:4)) %>%
  rename("Sample_Identifier" = "X.Sample.Identifier",
         "Patient_Identifier" = "Patient.Identifier",
         # "Other_Sample_ID" = "Other.Sample.ID",
         "Neuroendocrine_Features" = "Neuroendocrine.Features",
         "NEPC_Score" = "NEPC.Score",
         "Tissue_Site" = "Tissue.Site",
         "Pathology_Classification" = "Pathology.Classification",
         "Taxane_Exposure_Status" = "Taxane.Exposure.Status") %>%
  dplyr::select(Patient_Identifier,
                Sample_Identifier,
                # Other_Sample_ID,
                Neuroendocrine_Features,
                NEPC_Score,
                Tissue_Site,
                Pathology_Classification,
                Taxane_Exposure_Status)

# print("Combining the clinical sample and clinical patient data...")
data_combined <- left_join(data_clinical_sample,
                           data_clinical_patient,
                           by = c("Patient_Identifier" = "Patient_Identifier")) %>%
  filter(Sample_Identifier %in% colnames(gExpr)) # Filter metadata by that in the matrix

data_combined$Pathology_Classification <- data_combined$Pathology_Classification %>%
  stringr::str_replace_all("Not available", "Unclassified") %>%
  stringr::str_replace_all("Inadequate for diagnosis", "Unclassified") %>%
  stringr::str_replace_all("Small cell", "Pure small cell")
data_combined$NEPC_Score <- as.numeric(data_combined$NEPC_Score)
dir.create("experiments/survival-analysis/processed_data/metadata/", showWarnings = FALSE)
saveRDS(data_combined, "experiments/survival-analysis/processed_data/metadata/SU2C_EastCoast_metadata_forSA.rds")

su2cEC_metadata_toPlot <- data_combined %>%
  dplyr::select("Sample_Identifier",
                "Tissue_Site",
                "Pathology_Classification",
                "Neuroendocrine_Features",
                "NEPC_Score"
  )

saveRDS(su2cEC_metadata_toPlot, "experiments/oncomatch-analysis/processed_data/metadata/SU2C_EastCoast_metadata_forOMPlot.rds")
