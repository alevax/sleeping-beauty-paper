#!/usr/bin/env Rscript
#
# fusion_finder_analysis.R
#
# Consolidated fusion finder analysis pipeline:
#   1. CIS-Fusion overlap analysis
#   2. Figure S4B generation (box plot)
#   3. Figure S4B-C generation (box + bar plots)
#
# This script combines:
#   - overlap_cis_fusions.R
#   - figure_S4B_fusions.R
#   - figure_S4BC_fusions.R
#
# Prerequisites:
#   - HISAT2 alignment complete (results/bams/*.bam)
#   - FUSION_FINDER.pl complete (results/fusions/fusions_all.txt)
#
# Usage: Rscript libs/fusion-finder/fusion_finder_analysis.R [project_dir]
#

# === Load Libraries ===
suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(ggpubr)
  library(patchwork)
})

# === Configuration ===
args <- commandArgs(trailingOnly = TRUE)

if (length(args) > 0) {
  data_dir <- args[1]
} else {
  # Default to fusion_finder directory for input data
  data_dir <- "fusion_finder"
}

# Input paths (relative to repo root)
all_insertions_file <- "data/cis-results/all_insertions_master.csv"
key_genes_file <- "data/cis-results/key_genes_insertions.csv"
fusions_file <- file.path(data_dir, "results/fusions/fusions_all.txt")
samples_file <- file.path(data_dir, "data/fastq/samples.txt")
metadata_file <- file.path(data_dir, "R/RNAseq and CIS list Sep2020.xlsx")
bam_dir <- file.path(data_dir, "results/bams")

# Output directories (write to experiments/)
output_dir <- "experiments/fusion-finder"
figures_dir <- file.path(output_dir, "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

cat("=======================================================================\n")
cat("FUSION FINDER ANALYSIS PIPELINE\n")
cat("=======================================================================\n")
cat("Data directory:", data_dir, "\n")
cat("Output directory:", output_dir, "\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# =======================================================================
# PART 1: CIS-FUSION OVERLAP ANALYSIS
# =======================================================================

cat("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n")
cat("PART 1: CIS-Fusion Overlap Analysis\n")
cat("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n\n")

# === Load paired-end samples (fusion-capable) ===
paired_samples <- read_csv(samples_file, col_names = c("sample_id", "seq_type"),
                           show_col_types = FALSE) %>%
  filter(seq_type == "paired") %>%
  pull(sample_id)

cat("Paired-end samples:", length(paired_samples), "\n")

# === Load fusion data ===
if (!file.exists(fusions_file)) {
  stop("No fusions file found at ", fusions_file)
}

fusions <- read_tsv(fusions_file, col_names = c("chr", "start", "end", "sample",
                                                 "count1", "count2", "score", "chrSB_region"),
                    show_col_types = FALSE)
cat("Fusion events:", nrow(fusions), "\n")

# === Load CIS insertion data ===
all_insertions <- read_csv(all_insertions_file, show_col_types = FALSE)
key_genes <- read_csv(key_genes_file, show_col_types = FALSE)

cat("All insertions:", nrow(all_insertions), "\n")
cat("Key gene insertions:", nrow(key_genes), "\n\n")

# === Helper: parse comma-separated rnaseq_id ===
parse_rnaseq_ids <- function(rnaseq_id) {
  if (is.na(rnaseq_id)) return(character(0))
  str_split(rnaseq_id, ",")[[1]] %>% str_trim()
}

# === Expand and filter insertions to paired-end samples ===
prep_insertions <- function(insertions) {
  expanded <- insertions %>%
    rowwise() %>%
    mutate(rnaseq_ids_list = list(parse_rnaseq_ids(rnaseq_id))) %>%
    ungroup() %>%
    unnest(rnaseq_ids_list, keep_empty = TRUE)

  names(expanded)[names(expanded) == "rnaseq_ids_list"] <- "rnaseq_sample"

  expanded %>% filter(rnaseq_sample %in% paired_samples)
}

all_ins_pe <- prep_insertions(all_insertions)
key_ins_pe <- prep_insertions(key_genes)

# === Find overlaps ===
find_overlaps <- function(insertions, fusions) {
  if (nrow(insertions) == 0 || nrow(fusions) == 0) return(tibble())

  results <- list()

  for (i in seq_len(nrow(fusions))) {
    f <- fusions[i, ]

    matches <- insertions %>%
      filter(rnaseq_sample == f$sample) %>%
      mutate(
        overlaps_cis = (gene_chromosome == f$chr & cis_start <= f$end & cis_end >= f$start),
        overlaps_insertion = (insertion_chr == f$chr & insertion_start <= f$end & insertion_end >= f$start),
        overlaps_gene = (gene_chromosome == f$chr & gene_start <= f$end & gene_end >= f$start)
      ) %>%
      filter(overlaps_cis | overlaps_insertion | overlaps_gene) %>%
      mutate(
        overlap_type = case_when(
          overlaps_insertion ~ "insertion",
          overlaps_gene ~ "gene",
          overlaps_cis ~ "cis"
        ),
        fusion_chr = f$chr,
        fusion_start = f$start,
        fusion_end = f$end,
        fusion_score = f$score,
        fusion_chrSB_region = f$chrSB_region
      ) %>%
      select(-overlaps_cis, -overlaps_insertion, -overlaps_gene)

    if (nrow(matches) > 0) results[[length(results) + 1]] <- matches
  }

  if (length(results) == 0) return(tibble())
  bind_rows(results)
}

# === Run overlaps ===
cat("Finding overlaps...\n")
overlaps_all <- find_overlaps(all_ins_pe, fusions)
overlaps_key <- find_overlaps(key_ins_pe, fusions)

cat("  All genes overlaps:", nrow(overlaps_all), "\n")
cat("  Key genes overlaps:", nrow(overlaps_key), "\n\n")

# === Create clean output tables ===
format_output <- function(overlaps) {
  if (nrow(overlaps) == 0) return(tibble())

  overlaps %>%
    select(
      gene_name,
      sample = rnaseq_sample,
      is_nepc,
      overlap_type,
      insertion_chr,
      insertion_start,
      insertion_end,
      cis_start,
      cis_end,
      fusion_chr,
      fusion_start,
      fusion_end,
      fusion_score,
      fusion_chrSB_region,
      classification,
      sb_mechanism
    ) %>%
    distinct() %>%
    arrange(gene_name, sample)
}

out_all <- format_output(overlaps_all)
out_key <- format_output(overlaps_key)

# === Write overlap outputs ===
write_csv(out_all, file.path(output_dir, "all_genes_cis_fusions.csv"))
cat("Written: all_genes_cis_fusions.csv\n")

write_csv(out_key, file.path(output_dir, "key_genes_cis_fusions.csv"))
cat("Written: key_genes_cis_fusions.csv\n\n")

# === Sample Summary with Model Info ===
cat("Generating sample summary...\n")

# Load sample metadata
metadata <- read_excel(metadata_file, sheet = "RNAseq Min GA") %>%
  filter(!is.na(`RNAseq ID`)) %>%
  select(Model, `NEPC phenotype`, `RNAseq ID`) %>%
  rename(sample = `RNAseq ID`, nepc_phenotype = `NEPC phenotype`)

# Count fusions per sample
fusion_counts <- fusions %>%
  group_by(sample) %>%
  summarise(n_fusions = n(), .groups = "drop")

# Count CIS insertions per sample
cis_counts <- all_ins_pe %>%
  group_by(rnaseq_sample) %>%
  summarise(n_cis_insertions = n_distinct(paste(insertion_chr, insertion_start)), .groups = "drop") %>%
  rename(sample = rnaseq_sample)

# Get chrSB read counts from BAM files
cat("Counting chrSB reads from BAM files...\n")
chrSB_counts <- tibble(sample = paired_samples) %>%
  rowwise() %>%
  mutate(
    bam_file = file.path(bam_dir, paste0(sample, "_sorted.bam")),
    chrSB_reads = {
      if (file.exists(bam_file)) {
        cmd <- paste("samtools view -c", bam_file, "chrSB 2>/dev/null")
        result <- suppressWarnings(system(cmd, intern = TRUE))
        as.integer(result[1])
      } else {
        NA_integer_
      }
    }
  ) %>%
  ungroup() %>%
  select(sample, chrSB_reads)

# Build sample summary
sample_summary <- tibble(sample = paired_samples) %>%
  left_join(metadata, by = "sample") %>%
  left_join(fusion_counts, by = "sample") %>%
  left_join(cis_counts, by = "sample") %>%
  left_join(chrSB_counts, by = "sample") %>%
  mutate(
    n_fusions = replace_na(n_fusions, 0),
    n_cis_insertions = replace_na(n_cis_insertions, 0),
    chrSB_reads = replace_na(chrSB_reads, 0),
    has_transposon = str_detect(Model, "SBT2\\(\\+\\)"),
    nepc_phenotype = as.logical(nepc_phenotype)
  ) %>%
  arrange(Model, sample)

# Write sample summary
write_csv(sample_summary, file.path(output_dir, "sample_summary.csv"))
cat("Written: sample_summary.csv\n\n")

# === Generate comprehensive fusion tables ===
if (nrow(out_all) > 0) {
  cis_overlap_lookup <- out_all %>%
    select(sample, fusion_chr, fusion_start, fusion_end, gene_name, overlap_type, classification) %>%
    group_by(sample, fusion_chr, fusion_start, fusion_end) %>%
    summarise(
      cis_genes = paste(unique(gene_name), collapse = "; "),
      overlap_types = paste(unique(overlap_type), collapse = "; "),
      classifications = paste(unique(classification), collapse = "; "),
      .groups = "drop"
    )
} else {
  cis_overlap_lookup <- tibble()
}

fusion_table <- fusions %>%
  left_join(metadata, by = "sample") %>%
  {
    if (nrow(cis_overlap_lookup) > 0) {
      left_join(., cis_overlap_lookup, by = c("sample", "chr" = "fusion_chr", "start" = "fusion_start", "end" = "fusion_end"))
    } else {
      mutate(., cis_genes = NA_character_, overlap_types = NA_character_, classifications = NA_character_)
    }
  } %>%
  mutate(
    overlaps_cis = !is.na(cis_genes),
    chr_display = str_replace(chr, "^chr", ""),
    genomic_location = paste0("chr", chr_display, ":", start, "-", end),
    model_simple = case_when(
      str_detect(Model, "NPSBT2\\(\\+\\)") ~ "NP-SB(+)",
      str_detect(Model, "NPSBT2\\(-\\)") ~ "NP-SB(-)",
      str_detect(Model, "NPp53SBT2\\(\\+\\)") ~ "NPp53-SB(+)",
      str_detect(Model, "NPp53SBT2\\(-\\)") ~ "NPp53-SB(-)",
      TRUE ~ Model
    ),
    nepc_status = ifelse(nepc_phenotype == TRUE, "NEPC", "non-NEPC")
  ) %>%
  select(
    sample,
    model = model_simple,
    nepc_status,
    chr = chr_display,
    start,
    end,
    genomic_location,
    fusion_score = score,
    transposon_region = chrSB_region,
    read_count_1 = count1,
    read_count_2 = count2,
    overlaps_cis,
    cis_genes,
    overlap_type = overlap_types,
    classification = classifications
  ) %>%
  arrange(model, desc(nepc_status), chr, start)

write_csv(fusion_table, file.path(output_dir, "all_fusion_events_table.csv"))
cat("Written: all_fusion_events_table.csv\n")

if (sum(fusion_table$overlaps_cis) > 0) {
  cis_fusion_table <- fusion_table %>%
    filter(overlaps_cis) %>%
    arrange(cis_genes, sample)

  write_csv(cis_fusion_table, file.path(output_dir, "cis_overlapping_fusion_events_table.csv"))
  cat("Written: cis_overlapping_fusion_events_table.csv\n")
}

fusion_sample_summary <- fusion_table %>%
  group_by(sample, model, nepc_status) %>%
  summarise(
    n_total_fusions = n(),
    n_cis_fusions = sum(overlaps_cis),
    fusion_score_mean = round(mean(fusion_score), 1),
    fusion_score_median = median(fusion_score),
    transposon_regions = paste(sort(unique(transposon_region)), collapse = "; "),
    cis_genes = paste(unique(na.omit(cis_genes)), collapse = "; "),
    .groups = "drop"
  ) %>%
  mutate(cis_genes = ifelse(cis_genes == "", NA_character_, cis_genes)) %>%
  arrange(model, desc(nepc_status), desc(n_total_fusions))

write_csv(fusion_sample_summary, file.path(output_dir, "fusion_events_by_sample.csv"))
cat("Written: fusion_events_by_sample.csv\n\n")

cat("Part 1 complete: CIS-Fusion overlap analysis done\n\n")

# =======================================================================
# PART 2: FIGURE S4B (BOX PLOT)
# =======================================================================

cat("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n")
cat("PART 2: Figure S4B - Fusion Events by Tumor Type\n")
cat("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n\n")

# Prepare data for plotting
plot_data_s4b <- sample_summary %>%
  filter(str_detect(Model, "^NPp53")) %>%
  mutate(
    group = case_when(
      str_detect(Model, "SBT2\\(-\\)") ~ "NPp53-SB(-)",
      str_detect(Model, "SBT2\\(\\+\\)") & nepc_phenotype == TRUE ~ "NPp53-SB(+)\nNEPC",
      str_detect(Model, "SBT2\\(\\+\\)") & nepc_phenotype == FALSE ~ "NPp53-SB(+)\nnon-NEPC",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(group))

plot_data_s4b$group <- factor(plot_data_s4b$group,
                              levels = c("NPp53-SB(-)", "NPp53-SB(+)\nnon-NEPC", "NPp53-SB(+)\nNEPC"))

# Comparisons
comparisons_s4b <- list(
  c("NPp53-SB(-)", "NPp53-SB(+)\nnon-NEPC"),
  c("NPp53-SB(-)", "NPp53-SB(+)\nNEPC")
)

# Colors
colors_s4b <- c("NPp53-SB(-)" = "#808080",
                "NPp53-SB(+)\nnon-NEPC" = "#377EB8",
                "NPp53-SB(+)\nNEPC" = "#E41A1C")

# Create plot
p_s4b <- ggplot(plot_data_s4b, aes(x = group, y = n_fusions, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6) +
  geom_jitter(width = 0.15, size = 2.5, alpha = 0.8) +
  scale_fill_manual(values = colors_s4b) +
  stat_compare_means(comparisons = comparisons_s4b,
                     method = "wilcox.test",
                     label = "p.format",
                     tip.length = 0.02,
                     step.increase = 0.12,
                     size = 3.5) +
  labs(x = NULL, y = "Number of fusion events", title = NULL) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 11)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.20)))

# Save Figure (box plot)
ggsave(file.path(figures_dir, "fusions_by_tumor_type.pdf"), p_s4b, width = 4, height = 4, units = "in")
ggsave(file.path(figures_dir, "fusions_by_tumor_type.png"), p_s4b, width = 4, height = 4, units = "in", dpi = 300)
cat("Saved: fusions_by_tumor_type.pdf\n")
cat("Saved: fusions_by_tumor_type.png\n\n")

cat("Part 2 complete: Figure S4B generated\n\n")

# =======================================================================
# PART 3: FIGURE S4B-C (BOX + BAR PLOTS)
# =======================================================================

cat("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n")
cat("PART 3: Figure S4B-C - Combined Panels\n")
cat("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n\n")

# Panel A: Same as S4B (box plot by tumor type)
panel_A <- p_s4b

# Panel B: Fusion counts by transposon region
region_order <- c("5-flank", "SA1", "MSCV-LTR", "SD", "SA2", "IRDR", "3-flank")

sample_nepc <- sample_summary %>%
  select(sample, nepc_phenotype) %>%
  mutate(tumor_type = ifelse(nepc_phenotype, "NEPC", "non-NEPC"))

plot_data_s4c <- fusions %>%
  left_join(sample_nepc, by = "sample") %>%
  filter(!is.na(tumor_type)) %>%
  count(chrSB_region, tumor_type, name = "count") %>%
  rename(region = chrSB_region) %>%
  complete(region = region_order, tumor_type = c("non-NEPC", "NEPC"), fill = list(count = 0)) %>%
  mutate(
    region = factor(region, levels = region_order),
    tumor_type = factor(tumor_type, levels = c("non-NEPC", "NEPC"))
  )

colors_s4c <- c("non-NEPC" = "#377EB8", "NEPC" = "#E41A1C")

panel_B <- ggplot(plot_data_s4c, aes(x = region, y = count, fill = tumor_type)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = colors_s4c, name = "Tumor type") +
  labs(x = NULL, y = "Number of fusion events") +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    plot.margin = margin(5, 5, 5, 10)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))

# Save transposon region panel
ggsave(file.path(figures_dir, "fusions_by_transposon_region.pdf"), panel_B, width = 4, height = 3.5, units = "in")
ggsave(file.path(figures_dir, "fusions_by_transposon_region.png"), panel_B, width = 4, height = 3.5, units = "in", dpi = 300)
cat("Saved: fusions_by_transposon_region.pdf\n")
cat("Saved: fusions_by_transposon_region.png\n\n")

cat("Part 3 complete: Figure S4B-C panels generated\n\n")

# =======================================================================
# SUMMARY
# =======================================================================

cat("=======================================================================\n")
cat("FUSION FINDER ANALYSIS COMPLETE\n")
cat("=======================================================================\n")
cat("End time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("Output files:\n")
cat("  Overlap tables:\n")
cat("    - all_genes_cis_fusions.csv\n")
cat("    - key_genes_cis_fusions.csv\n")
cat("    - sample_summary.csv\n")
cat("    - all_fusion_events_table.csv\n")
cat("    - fusion_events_by_sample.csv\n")
cat("  Figures:\n")
cat("    - fusions_by_tumor_type.pdf/png\n")
cat("    - fusions_by_transposon_region.pdf/png\n")
cat("\n=======================================================================\n")
