# SIRT1 Expression and Activity Analysis
# Compares NEPC vs Adeno in Beltran 2016 and SU2C EastCoast datasets

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(ggpubr)
  library(viper)
})

output_dir <- "experiments/sirt1-human-expression"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Define color palette
group_colors <- c("Adeno" = "#4DBBD5", "NEPC" = "#E64B35")

# ==============================================================================
# BELTRAN 2016 DATASET
# ==============================================================================
cat("=== SIRT1 Analysis in Beltran 2016 Dataset ===\n\n")

beltran <- readRDS("data/Beltran-2016/data_obj_Beltran_non-NEPC_vs_NEPC.rds")
metadata_beltran <- beltran$metadata
expr_mat_beltran <- beltran$data
viper_mat_beltran <- beltran$vpmat

# Get SIRT1 data
sirt1_expr_beltran <- unlist(expr_mat_beltran["SIRT1", ])
names(sirt1_expr_beltran) <- colnames(expr_mat_beltran)
sirt1_viper_beltran <- unlist(viper_mat_beltran["SIRT1", ])
names(sirt1_viper_beltran) <- colnames(viper_mat_beltran)

# Build plotting data frame
plot_df_beltran <- metadata_beltran %>%
  mutate(
    SIRT1_expression = sirt1_expr_beltran[sample_id],
    SIRT1_activity = sirt1_viper_beltran[sample_id],
    Group = case_when(
      histological_type == "CRPC-NE" ~ "NEPC",
      histological_type == "CRPC-Adeno" ~ "Adeno",
      TRUE ~ "Other"
    ),
    Dataset = "Beltran 2016"
  ) %>%
  filter(Group %in% c("NEPC", "Adeno"))

plot_df_beltran$Group <- factor(plot_df_beltran$Group, levels = c("Adeno", "NEPC"))

cat("Sample counts:\n")
print(table(plot_df_beltran$Group))

expr_test_beltran <- t.test(SIRT1_expression ~ Group, data = plot_df_beltran)
viper_test_beltran <- t.test(SIRT1_activity ~ Group, data = plot_df_beltran)

cat("\n--- Expression Analysis ---\n")
cat("Mean SIRT1 expression:\n")
print(tapply(plot_df_beltran$SIRT1_expression, plot_df_beltran$Group, mean, na.rm = TRUE))
cat(sprintf("Welch's t-test p-value: %.4f\n", expr_test_beltran$p.value))

cat("\n--- VIPER Activity Analysis ---\n")
cat("Mean SIRT1 activity:\n")
print(tapply(plot_df_beltran$SIRT1_activity, plot_df_beltran$Group, mean, na.rm = TRUE))
cat(sprintf("Welch's t-test p-value: %.4f\n", viper_test_beltran$p.value))

# ==============================================================================
# SU2C EASTCOAST DATASET
# ==============================================================================
cat("\n\n=== SIRT1 Analysis in SU2C EastCoast Dataset ===\n\n")

# Load the processed RDS file which has harmonized sample IDs
su2c_data <- readRDS("data/SU2C_EastCoast/abida-2019-processed-data.rds")
su2c_expr <- su2c_data$counts_fpkm.polyA
su2c_meta <- su2c_data$metadata

# Handle duplicate gene names by keeping first occurrence
su2c_expr <- su2c_expr[!duplicated(su2c_expr$Hugo_Symbol), ]
rownames(su2c_expr) <- su2c_expr$Hugo_Symbol
su2c_expr <- su2c_expr[, -1]  # Remove Hugo_Symbol column

# Load interactome for VIPER
su2c_net <- readRDS("data/SU2C_EastCoast/SU2CEastCoast_net_pruned.rds")

# Get SIRT1 expression
sirt1_expr_su2c <- as.numeric(su2c_expr["SIRT1", ])
names(sirt1_expr_su2c) <- colnames(su2c_expr)

# Run VIPER on SU2C data
cat("Running VIPER on SU2C EastCoast data...\n")
su2c_expr_mat <- as.matrix(su2c_expr)
su2c_viper <- viper(su2c_expr_mat, su2c_net, verbose = FALSE)

# Get SIRT1 VIPER activity
sirt1_viper_su2c <- su2c_viper["SIRT1", ]

# Get expression sample IDs
expr_samples <- colnames(su2c_expr)

# Filter metadata to samples with expression data
su2c_meta_matched <- su2c_meta[su2c_meta$SAMPLE_ID %in% expr_samples, ]

# Build plotting data frame
plot_df_su2c <- data.frame(
  sample_id = expr_samples,
  SIRT1_expression = sirt1_expr_su2c,
  SIRT1_activity = sirt1_viper_su2c[expr_samples],
  stringsAsFactors = FALSE
)

# Match to clinical data using SAMPLE_ID directly (already harmonized in RDS)
plot_df_su2c <- plot_df_su2c %>%
  left_join(su2c_meta_matched %>% select(SAMPLE_ID, NEUROENDOCRINE_FEATURES),
            by = c("sample_id" = "SAMPLE_ID"))

# Classify groups
plot_df_su2c <- plot_df_su2c %>%
  mutate(
    Group = case_when(
      NEUROENDOCRINE_FEATURES == "Yes" ~ "NEPC",
      NEUROENDOCRINE_FEATURES == "No" ~ "Adeno",
      TRUE ~ "Other"
    ),
    Dataset = "SU2C EastCoast"
  ) %>%
  filter(Group %in% c("NEPC", "Adeno"))

plot_df_su2c$Group <- factor(plot_df_su2c$Group, levels = c("Adeno", "NEPC"))

cat("Sample counts:\n")
print(table(plot_df_su2c$Group))

expr_test_su2c <- t.test(SIRT1_expression ~ Group, data = plot_df_su2c)
viper_test_su2c <- t.test(SIRT1_activity ~ Group, data = plot_df_su2c)

cat("\n--- Expression Analysis ---\n")
cat("Mean SIRT1 expression:\n")
print(tapply(plot_df_su2c$SIRT1_expression, plot_df_su2c$Group, mean, na.rm = TRUE))
cat(sprintf("Welch's t-test p-value: %.4g\n", expr_test_su2c$p.value))

cat("\n--- VIPER Activity Analysis ---\n")
cat("Mean SIRT1 activity:\n")
print(tapply(plot_df_su2c$SIRT1_activity, plot_df_su2c$Group, mean, na.rm = TRUE))
cat(sprintf("Welch's t-test p-value: %.4g\n", viper_test_su2c$p.value))

# ==============================================================================
# GENERATE PLOTS
# ==============================================================================

# Beltran plots
p_expr_beltran <- ggplot(plot_df_beltran, aes(x = Group, y = SIRT1_expression, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 1.5, alpha = 0.6) +
  scale_fill_manual(values = group_colors) +
  stat_compare_means(method = "t.test", label = "p.format",
                     label.x.npc = "center", label.y.npc = "top") +
  labs(title = "SIRT1 Expression", subtitle = "Beltran 2016 (n=49)", x = "", y = "SIRT1 Expression (CPM)") +
  theme_classic() +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))

p_viper_beltran <- ggplot(plot_df_beltran, aes(x = Group, y = SIRT1_activity, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 1.5, alpha = 0.6) +
  scale_fill_manual(values = group_colors) +
  stat_compare_means(method = "t.test", label = "p.format",
                     label.x.npc = "center", label.y.npc = "top") +
  labs(title = "SIRT1 Activity (VIPER)", subtitle = "Beltran 2016 (n=49)", x = "", y = "SIRT1 Activity (NES)") +
  theme_classic() +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))

# SU2C plots
p_expr_su2c <- ggplot(plot_df_su2c, aes(x = Group, y = SIRT1_expression, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 1.5, alpha = 0.6) +
  scale_fill_manual(values = group_colors) +
  stat_compare_means(method = "t.test", label = "p.format",
                     label.x.npc = "center", label.y.npc = "top") +
  labs(title = "SIRT1 Expression", subtitle = paste0("SU2C EastCoast (n=", nrow(plot_df_su2c), ")"),
       x = "", y = "SIRT1 Expression (FPKM)") +
  theme_classic() +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))

p_viper_su2c <- ggplot(plot_df_su2c, aes(x = Group, y = SIRT1_activity, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 1.5, alpha = 0.6) +
  scale_fill_manual(values = group_colors) +
  stat_compare_means(method = "t.test", label = "p.format",
                     label.x.npc = "center", label.y.npc = "top") +
  labs(title = "SIRT1 Activity (VIPER)", subtitle = paste0("SU2C EastCoast (n=", nrow(plot_df_su2c), ")"),
       x = "", y = "SIRT1 Activity (NES)") +
  theme_classic() +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))

# Save individual plots only (no combined plot with A/B/C/D labels)
ggsave(file.path(output_dir, "SIRT1_expression_beltran.pdf"), p_expr_beltran, width = 4, height = 5)
ggsave(file.path(output_dir, "SIRT1_viper_beltran.pdf"), p_viper_beltran, width = 4, height = 5)
ggsave(file.path(output_dir, "SIRT1_expression_su2c.pdf"), p_expr_su2c, width = 4, height = 5)
ggsave(file.path(output_dir, "SIRT1_viper_su2c.pdf"), p_viper_su2c, width = 4, height = 5)

# Save source data
write.csv(plot_df_beltran, file.path(output_dir, "SIRT1_beltran_source_data.csv"), row.names = FALSE)
write.csv(plot_df_su2c, file.path(output_dir, "SIRT1_su2c_source_data.csv"), row.names = FALSE)

cat("\n=== Output files saved to", output_dir, "===\n")
cat("- SIRT1_expression_beltran.pdf\n")
cat("- SIRT1_viper_beltran.pdf\n")
cat("- SIRT1_expression_su2c.pdf\n")
cat("- SIRT1_viper_su2c.pdf\n")
cat("- SIRT1_beltran_source_data.csv\n")
cat("- SIRT1_su2c_source_data.csv\n")
