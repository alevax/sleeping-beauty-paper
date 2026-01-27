#!/usr/bin/env Rscript
#
# chimeric_isoforms_analysis.R
#
# Consolidated analysis pipeline for chimeric isoform detection in Sleeping Beauty
# transposon insertions from RNA-seq data. This script combines three analysis stages:
#
#   Part 1: Data Preparation - Filters CIS insertions for NEPC samples with RNA-seq
#   Part 2: Coverage Analysis - Extracts RNA-seq coverage and computes ratios
#   Part 3: Publication Plots - Generates publication-quality coverage visualizations
#
# Input requirements:
#   - data/reference/key_genes_insertions.csv (CIS results)
#   - data/reference/RNAseq and CIS list Sep2020.xlsx (sample metadata)
#   - data/reference/gencode.vM25.annotation.gtf (gene annotations)
#   - data/bams/*.bam (aligned RNA-seq BAM files)
#
# Output files:
#   - data/chimeric_data.rds (intermediate data object)
#   - data/chimeric_genes_insertions.csv (filtered insertions)
#   - data/sample_groups.csv (sample classifications)
#   - results/coverage_data/{gene}_coverage.csv (per-gene coverage data)
#   - results/plots/{gene}_coverage.pdf (raw coverage plots)
#   - results/plots/{gene}_coverage_normalized.pdf (normalized coverage plots)
#   - results/tables/coverage_ratios.csv (upstream/downstream ratios)
#   - results/tables/analysis_summary.csv (gene-level summary)
#
# Usage: Rscript chimeric_isoforms_analysis.R
#
# Author: Matteo DiBernardo
# Date: 2026-01-26
#

# =============================================================================
# LIBRARY LOADING
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)        # Data manipulation and visualization
  library(readxl)          # Excel file reading
  library(rtracklayer)     # GTF/GFF file handling
  library(Rsamtools)       # BAM file operations
  library(GenomicAlignments) # RNA-seq alignment processing
  library(GenomicFeatures) # Gene structure annotations
  library(scales)          # Plot formatting
  library(patchwork)       # Multi-panel plot layouts
})

# =============================================================================
# GLOBAL CONFIGURATION
# =============================================================================

# Input files (relative to repo root)
insertions_file <- "data/cis-results/key_genes_insertions.csv"
sample_xlsx <- "data/gene-expression-from-sb-mouse-models/RNAseq and CIS list Sep2020.xlsx"
gtf_file <- "isoform_switching/reference/gencode.vM25.annotation.gtf"
bam_dir <- "fusion_finder/results/bams"

# Output directories (write to experiments/)
output_dir <- "experiments/chimeric-isoforms"
coverage_dir <- file.path(output_dir, "coverage_data")
plots_dir <- file.path(output_dir, "plots")
tables_dir <- file.path(output_dir, "tables")
processed_data_dir <- file.path(output_dir, "processed_data")

# Output files
output_insertions <- file.path(processed_data_dir, "chimeric_genes_insertions.csv")
output_samples <- file.path(processed_data_dir, "sample_groups.csv")
output_rds <- file.path(processed_data_dir, "chimeric_data.rds")

# Analysis parameters
BIN_SIZE <- 100
FLANK_SIZE <- 10000

# Color palettes
COLORS <- list(
  nepc = "#C73E1D",           # Red-brown for NEPC
  sb_neg = "#2E86AB",         # Steel blue for SB-negative
  sb_nonne = "#F18F01",       # Orange for SB-positive nonNE
  exon = "#2c3e50",           # Dark slate for exons
  intron = "#95a5a6",         # Gray for intron line
  sense = "#E74C3C",          # Red for sense orientation
  antisense = "#3498DB",      # Blue for antisense orientation
  mixed = "#9B59B6"           # Purple for mixed orientation
)

cat("=======================================================================\n")
cat("CHIMERIC ISOFORMS ANALYSIS PIPELINE\n")
cat("=======================================================================\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("BAM directory:", bam_dir, "\n")
cat("Output directory:", output_dir, "\n\n")

# =============================================================================
# PART 1: DATA PREPARATION
# =============================================================================

cat("\n")
cat("=======================================================================\n")
cat("PART 1: DATA PREPARATION\n")
cat("=======================================================================\n\n")

# Create output directories for Part 1
dir.create(processed_data_dir, recursive = TRUE, showWarnings = FALSE)

cat("Filtering CIS insertions for:\n")
cat("  - NEPC samples (is_nepc == TRUE)\n")
cat("  - Has RNA-seq data (has_rnaseq == TRUE)\n")
cat("  - BAM file exists\n")
cat("  - Insertion within gene body (verified from GTF)\n\n")

# === Step 1: Get Available BAM Samples ===
cat("Step 1: Identifying available BAM files...\n")

bam_files <- list.files(bam_dir, pattern = "_sorted\\.bam$", full.names = FALSE)
available_samples <- str_remove(bam_files, "_sorted\\.bam")

cat("  BAM files found:", length(available_samples), "\n")
cat("  Samples:", paste(head(available_samples, 5), collapse = ", "), "...\n\n")

# === Step 2: Load GTF for Gene Boundaries ===
cat("Step 2: Loading GTF for gene boundaries...\n")

gtf <- rtracklayer::import(gtf_file)
genes_gtf <- gtf[gtf$type == "gene"]

# Create gene boundary lookup
gene_boundaries <- as.data.frame(genes_gtf) %>%
  dplyr::select(gene_name, seqnames, start, end, strand) %>%
  dplyr::rename(
    gtf_chr = seqnames,
    gtf_start = start,
    gtf_end = end,
    gtf_strand = strand
  ) %>%
  filter(!is.na(gene_name))

cat("  Genes in GTF:", nrow(gene_boundaries), "\n\n")

# === Step 3: Load and Parse Sample Metadata ===
cat("Step 3: Loading sample metadata from xlsx...\n")

rnaseq_meta <- read_excel(sample_xlsx, sheet = "RNAseq Min GA")

# Clean column names and remove corrupted first row
rnaseq_meta <- rnaseq_meta %>%
  rename_all(~ str_replace_all(., " ", "_")) %>%
  filter(!is.na(RNAseq_ID), Model != "43313")

cat("  Total RNA-seq samples in metadata:", nrow(rnaseq_meta), "\n")

# Classify samples by SB status
sample_groups <- rnaseq_meta %>%
  transmute(
    sample = RNAseq_ID,
    model = Model,
    nepc = as.logical(NEPC_phenotype),
    sb_negative = str_detect(Model, "SBT2\\(-\\)"),
    sb_positive = str_detect(Model, "SBT2\\(\\+\\)"),
    group = case_when(
      sb_negative ~ "SB_negative_control",
      sb_positive & nepc ~ "SB_positive_NEPC",
      sb_positive & !nepc ~ "SB_positive_nonNE",
      TRUE ~ "other"
    )
  ) %>%
  filter(sample %in% available_samples)

cat("  Samples with available BAMs:", nrow(sample_groups), "\n")
cat("\n  Sample groups:\n")
print(table(sample_groups$group))

write_csv(sample_groups, output_samples)
cat("\n  Saved:", output_samples, "\n\n")

# === Step 4: Load and Filter CIS Insertions ===
cat("Step 4: Loading and filtering CIS insertions...\n")

insertions_raw <- read_csv(insertions_file, show_col_types = FALSE)
cat("  Raw insertions:", nrow(insertions_raw), "\n")
cat("  Genes:", length(unique(insertions_raw$gene_name)), "\n")

# Filter 1 & 2: NEPC + has RNA-seq
insertions_nepc <- insertions_raw %>%
  filter(is_nepc == TRUE, has_rnaseq == TRUE)

cat("\n  After NEPC + has_rnaseq filter:", nrow(insertions_nepc), "\n")

# Expand rnaseq_id (may be comma-separated) and filter for available BAMs
insertions_expanded <- insertions_nepc %>%
  mutate(sample_id_list = str_split(rnaseq_id, ",")) %>%
  unnest(sample_id_list) %>%
  mutate(sample_id = str_trim(sample_id_list)) %>%
  filter(sample_id %in% available_samples) %>%
  dplyr::select(-sample_id_list)

cat("  After BAM availability filter:", nrow(insertions_expanded), "\n")

# Filter 4: Insertion within gene body (recompute from GTF)
insertions_with_gtf <- insertions_expanded %>%
  left_join(gene_boundaries, by = "gene_name")

# Check if insertion midpoint falls within gene boundaries
insertions_in_gene <- insertions_with_gtf %>%
  mutate(
    in_gene_body = insertion_midpoint >= gtf_start & insertion_midpoint <= gtf_end
  ) %>%
  filter(in_gene_body == TRUE)

cat("  After in-gene-body filter (from GTF):", nrow(insertions_in_gene), "\n")
cat("  Unique genes:", length(unique(insertions_in_gene$gene_name)), "\n")
cat("  Unique samples:", length(unique(insertions_in_gene$sample_id)), "\n")

# Show what was filtered out
filtered_out <- insertions_with_gtf %>%
  mutate(in_gene_body = insertion_midpoint >= gtf_start & insertion_midpoint <= gtf_end) %>%
  filter(in_gene_body == FALSE)

if (nrow(filtered_out) > 0) {
  cat("\n  Insertions filtered out (outside gene body):\n")
  print(filtered_out %>% dplyr::select(gene_name, sample_id, insertion_midpoint, gtf_start, gtf_end) %>% distinct())
}

# Summarize by gene
cat("\n  Gene summary:\n")
gene_summary <- insertions_in_gene %>%
  group_by(gene_name) %>%
  summarise(
    n_insertions = n(),
    n_samples = n_distinct(sample_id),
    samples = paste(unique(sample_id), collapse = ", "),
    .groups = "drop"
  ) %>%
  arrange(desc(n_insertions))

print(gene_summary, n = 20)

# === Step 5: Save Filtered Insertions ===
cat("\nStep 5: Saving filtered insertions...\n")

chimeric_insertions <- insertions_in_gene %>%
  dplyr::select(
    gene_name,
    sample_id,
    gene_chromosome,
    gene_strand,
    gene_start,
    gene_end,
    insertion_chr,
    insertion_start,
    insertion_end,
    insertion_midpoint,
    insertion_strand,
    classification,
    sb_mechanism,
    intron_number,
    pct_into_transcript,
    is_early_intronic,
    is_exonic,
    is_same_strand,
    gtf_start,
    gtf_end
  ) %>%
  left_join(
    sample_groups %>% dplyr::select(sample, group, nepc),
    by = c("sample_id" = "sample")
  ) %>%
  arrange(gene_name, sample_id)

write_csv(chimeric_insertions, output_insertions)
cat("  Saved:", output_insertions, "\n")

# === Step 6: Control Sample Summary ===
cat("\nStep 6: Preparing control sample groups...\n")

sb_neg_controls <- sample_groups %>%
  filter(group == "SB_negative_control") %>%
  pull(sample)

sb_pos_nonne_controls <- sample_groups %>%
  filter(group == "SB_positive_nonNE") %>%
  pull(sample)

# Save as RDS for faster loading in subsequent steps
saveRDS(
  list(
    insertions = chimeric_insertions,
    sample_groups = sample_groups,
    gene_boundaries = gene_boundaries,
    sb_neg_controls = sb_neg_controls,
    sb_pos_nonne_controls = sb_pos_nonne_controls,
    available_samples = available_samples
  ),
  output_rds
)
cat("  Saved RDS:", output_rds, "\n")

cat("  SB-negative controls (no transposon):", length(sb_neg_controls), "\n")
cat("    ", paste(sb_neg_controls, collapse = ", "), "\n")
cat("  SB-positive nonNE controls:", length(sb_pos_nonne_controls), "\n")
cat("    ", paste(sb_pos_nonne_controls, collapse = ", "), "\n")

cat("\nPart 1 Complete - Filtered dataset summary:\n")
cat("  Genes:", length(unique(chimeric_insertions$gene_name)), "\n")
cat("  Insertion-sample pairs:", nrow(chimeric_insertions), "\n")
cat("  Samples with insertions:", length(unique(chimeric_insertions$sample_id)), "\n")

# =============================================================================
# PART 2: COVERAGE ANALYSIS
# =============================================================================

cat("\n")
cat("=======================================================================\n")
cat("PART 2: COVERAGE ANALYSIS\n")
cat("=======================================================================\n\n")

cat("Extracting RNA-seq coverage and computing upstream/downstream ratios\n")
cat("Comparing NEPC samples with insertions to two control groups:\n")
cat("  1. SB-negative: NPp53SBT2(-) - no transposon\n")
cat("  2. SB-positive nonNE: NPp53SBT2(+)-nonNE - transposon but non-NEPC phenotype\n\n")

# === Create Output Directories ===
dir.create(processed_data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(coverage_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

# === Load Prepared Data from Part 1 ===
cat("Loading prepared data from Part 1...\n")

data <- readRDS(output_rds)
insertions <- data$insertions
sample_groups <- data$sample_groups
sb_neg_controls <- data$sb_neg_controls
sb_pos_nonne_controls <- data$sb_pos_nonne_controls

cat("  Genes:", length(unique(insertions$gene_name)), "\n")
cat("  Insertion-sample pairs:", nrow(insertions), "\n")
cat("  Samples with insertions:", length(unique(insertions$sample_id)), "\n")
cat("  SB-negative controls:", length(sb_neg_controls), "\n")
cat("  SB-positive nonNE controls:", length(sb_pos_nonne_controls), "\n\n")

# === Helper Functions ===

get_coverage <- function(bam_file, chr, start, end, bin_size = BIN_SIZE) {
  if (!file.exists(bam_file)) return(NULL)

  gr <- GRanges(seqnames = chr, ranges = IRanges(start = start, end = end))
  param <- ScanBamParam(which = gr)
  aln <- tryCatch(readGAlignments(bam_file, param = param), error = function(e) NULL)

  if (is.null(aln) || length(aln) == 0) {
    return(tibble(pos = seq(start, end, by = bin_size), coverage = 0))
  }

  cov <- coverage(aln)
  if (chr %in% names(cov)) {
    cov_vec <- as.numeric(cov[[chr]][start:min(end, length(cov[[chr]]))])
    if (length(cov_vec) < (end - start + 1)) {
      cov_vec <- c(cov_vec, rep(0, (end - start + 1) - length(cov_vec)))
    }
  } else {
    cov_vec <- rep(0, end - start + 1)
  }

  positions <- seq(start, end, by = bin_size)
  binned_cov <- sapply(positions, function(p) {
    idx_start <- p - start + 1
    idx_end <- min(idx_start + bin_size - 1, length(cov_vec))
    mean(cov_vec[idx_start:idx_end], na.rm = TRUE)
  })

  tibble(pos = positions, coverage = binned_cov)
}

calculate_ratio <- function(cov_df, insertion_pos, gene_strand, window = 5000) {
  if (gene_strand == "+") {
    up <- filter(cov_df, pos < insertion_pos, pos > insertion_pos - window)
    down <- filter(cov_df, pos > insertion_pos, pos < insertion_pos + window)
  } else {
    up <- filter(cov_df, pos > insertion_pos, pos < insertion_pos + window)
    down <- filter(cov_df, pos < insertion_pos, pos > insertion_pos - window)
  }

  up_mean <- mean(up$coverage, na.rm = TRUE)
  down_mean <- mean(down$coverage, na.rm = TRUE)

  if (is.na(up_mean) || up_mean == 0) {
    return(tibble(upstream = up_mean, downstream = down_mean, ratio = NA, log2_ratio = NA))
  }

  ratio <- down_mean / up_mean
  tibble(upstream = up_mean, downstream = down_mean, ratio = ratio,
         log2_ratio = if(ratio > 0) log2(ratio) else NA)
}

# === Process Each Gene ===
target_genes <- unique(insertions$gene_name)
cat("Processing", length(target_genes), "genes:\n")
cat("  ", paste(target_genes, collapse = ", "), "\n")
cat("-----------------------------------------------------------------------\n\n")

all_ratios <- tibble()
gene_summary_cov <- tibble()

for (gene_name in target_genes) {
  cat("=== ", gene_name, " ===\n", sep = "")

  # Gene info from insertions
  gene_data <- insertions %>% dplyr::filter(gene_name == !!gene_name) %>% dplyr::slice(1)
  chr <- gene_data$gene_chromosome
  strand <- gene_data$gene_strand
  gene_start <- gene_data$gtf_start
  gene_end <- gene_data$gtf_end

  plot_start <- max(1, gene_start - FLANK_SIZE)
  plot_end <- gene_end + FLANK_SIZE

  cat("  Location: ", chr, ":", comma(gene_start), "-", comma(gene_end), " (", strand, ")\n", sep = "")

  # Get exons from GTF
  exons <- gtf[gtf$gene_name == gene_name & gtf$type == "exon"]
  exons_df <- as.data.frame(exons) %>% dplyr::select(start, end) %>% distinct()
  cat("  Exons:", nrow(exons_df), "\n")

  # Get insertions for this gene
  gene_insertions <- insertions %>% filter(gene_name == !!gene_name)
  samples_with_insertion <- unique(gene_insertions$sample_id)
  cat("  Samples with insertion:", paste(samples_with_insertion, collapse = ", "), "\n")

  # Samples: insertion samples + both control groups
  samples_to_analyze <- unique(c(samples_with_insertion, sb_neg_controls, sb_pos_nonne_controls))
  cat("  Extracting coverage from", length(samples_to_analyze), "samples...\n")
  cat("    NEPC with insertion:", length(samples_with_insertion), "\n")
  cat("    SB-negative controls:", length(sb_neg_controls), "\n")
  cat("    SB-pos nonNE controls:", length(sb_pos_nonne_controls), "\n")

  # Extract coverage
  coverage_data <- map_dfr(samples_to_analyze, function(s) {
    bam_file <- file.path(bam_dir, paste0(s, "_sorted.bam"))
    cov <- get_coverage(bam_file, chr, plot_start, plot_end)
    if (!is.null(cov)) {
      cov %>%
        mutate(sample = s, has_insertion = s %in% samples_with_insertion) %>%
        left_join(sample_groups, by = "sample")
    }
  })

  if (nrow(coverage_data) == 0) {
    cat("  No coverage data, skipping.\n\n")
    next
  }

  write_csv(coverage_data, file.path(coverage_dir, paste0(gene_name, "_coverage.csv")))
  cat("  Saved coverage data: ", file.path(coverage_dir, paste0(gene_name, "_coverage.csv")), "\n")

  # Insertion coordinates
  insertion_coords <- gene_insertions %>%
    dplyr::select(insertion_start, insertion_end, insertion_midpoint) %>%
    distinct()

  ref_insertion_pos <- insertion_coords$insertion_midpoint[1]

  # === Calculate Ratios ===
  sample_ratios <- coverage_data %>%
    group_by(sample, has_insertion, group) %>%
    group_modify(~ calculate_ratio(.x, ref_insertion_pos, strand)) %>%
    ungroup() %>%
    mutate(gene_name = gene_name)

  all_ratios <- bind_rows(all_ratios, sample_ratios)

  # Calculate summary with both control groups
  mean_log2_nepc <- mean(sample_ratios$log2_ratio[sample_ratios$has_insertion], na.rm = TRUE)
  mean_log2_sbneg <- mean(sample_ratios$log2_ratio[sample_ratios$group == "SB_negative_control"], na.rm = TRUE)
  mean_log2_nonne <- mean(sample_ratios$log2_ratio[sample_ratios$group == "SB_positive_nonNE"], na.rm = TRUE)

  gene_summary_cov <- bind_rows(gene_summary_cov, tibble(
    gene_name = gene_name,
    chr = chr,
    strand = strand,
    n_samples_insertion = length(samples_with_insertion),
    n_sb_neg_controls = sum(sample_ratios$group == "SB_negative_control"),
    n_sb_nonne_controls = sum(sample_ratios$group == "SB_positive_nonNE"),
    mean_log2_NEPC = mean_log2_nepc,
    mean_log2_SBneg = mean_log2_sbneg,
    mean_log2_nonNE = mean_log2_nonne
  ))

  cat("  Mean log2 ratio (up/down):\n")
  cat("    NEPC with insertion:", round(mean_log2_nepc, 2), "\n")
  cat("    SB-negative control:", round(mean_log2_sbneg, 2), "\n")
  cat("    SB+ nonNE control:", round(mean_log2_nonne, 2), "\n")
  cat("  Done.\n\n")
}

# === Save Summary Tables ===
cat("-----------------------------------------------------------------------\n")
cat("Saving summary tables...\n")

write_csv(all_ratios, file.path(tables_dir, "coverage_ratios.csv"))
write_csv(gene_summary_cov, file.path(tables_dir, "analysis_summary.csv"))

cat("  Saved:", file.path(tables_dir, "coverage_ratios.csv"), "\n")
cat("  Saved:", file.path(tables_dir, "analysis_summary.csv"), "\n")

cat("\nPart 2 Complete - Results summary:\n")
if (nrow(gene_summary_cov) > 0 && "mean_log2_NEPC" %in% colnames(gene_summary_cov)) {
  print(gene_summary_cov %>% dplyr::select(gene_name, mean_log2_NEPC, mean_log2_SBneg, mean_log2_nonNE))
  cat("\nInterpretation:\n")
  cat("  Positive = higher downstream than upstream (normal transcription)\n")
  cat("  Negative = higher upstream (potential truncation)\n")
  cat("  Compare NEPC samples with insertion to both control groups:\n")
  cat("    SBneg = NPp53SBT2(-) - no transposon\n")
  cat("    nonNE = NPp53SBT2(+)-nonNE - transposon but non-NEPC\n")
} else {
  cat("  No coverage data available for analysis.\n")
}

cat("\nInterpretation (when data is available):\n")
cat("  Positive ratio = higher downstream than upstream coverage (normal transcription)\n")
cat("  Negative ratio = higher upstream coverage (potential truncation)\n")
cat("  SBneg = NPp53SBT2(-) control - no transposon\n")
cat("  nonNE = NPp53SBT2(+)-nonNE control - transposon but non-NEPC\n")

# =============================================================================
# PART 3: PUBLICATION PLOTS
# =============================================================================

cat("\n")
cat("=======================================================================\n")
cat("PART 3: PUBLICATION PLOTS\n")
cat("=======================================================================\n\n")

cat("Generating publication-quality coverage plots with:\n")
cat("  - Multi-panel layout (NEPC + controls)\n")
cat("  - Gene structure with exon/intron visualization\n")
cat("  - Lollipop markers for insertion sites\n")
cat("  - Both raw and normalized coverage versions\n\n")

# === Load GTF for Gene Structure ===
gtf_df <- as.data.frame(gtf)

#' Get exons for a gene from GTF
get_gene_exons <- function(gene_symbol, gtf_df) {
  exons <- gtf_df %>%
    dplyr::filter(type == "exon", gene_name == gene_symbol) %>%
    dplyr::select(seqnames, start, end, strand, transcript_id, exon_number) %>%
    dplyr::distinct()

  if (nrow(exons) == 0) return(NULL)

  # Collapse overlapping exons to get union
  exons_collapsed <- exons %>%
    dplyr::arrange(start) %>%
    dplyr::mutate(group = cumsum(c(1, diff(start) > 1))) %>%
    dplyr::group_by(seqnames, strand, group) %>%
    dplyr::summarise(start = min(start), end = max(end), .groups = "drop") %>%
    dplyr::select(-group)

  return(exons_collapsed)
}

#' Generate publication plot for a gene
generate_pub_plot <- function(gene_symbol) {

  cat("\n=== Generating publication plot for", gene_symbol, "===\n")

  # Load coverage data
  coverage_file <- file.path(coverage_dir, paste0(gene_symbol, "_coverage.csv"))
  if (!file.exists(coverage_file)) {
    cat("ERROR: Coverage file not found\n")
    return(NULL)
  }

  coverage_data <- read.csv(coverage_file)
  cat("Loaded", nrow(coverage_data), "coverage points\n")

  # Get gene info
  gene_info <- insertions %>%
    dplyr::filter(gene_name == gene_symbol) %>%
    dplyr::slice(1)

  chr <- gene_info$gene_chromosome
  gene_start <- gene_info$gene_start
  gene_end <- gene_info$gene_end
  gene_strand <- gene_info$gene_strand

  # Get insertions
  gene_insertions <- insertions %>%
    dplyr::filter(gene_name == gene_symbol)

  nepc_samples <- unique(gene_insertions$sample_id)
  cat("NEPC samples with insertions:", paste(nepc_samples, collapse = ", "), "\n")

  # === Compute averaged coverage ===

  nepc_coverage <- coverage_data %>%
    dplyr::filter(sample %in% nepc_samples) %>%
    dplyr::group_by(pos) %>%
    dplyr::summarise(coverage = mean(coverage, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(group = "NEPC")

  sb_neg_coverage <- coverage_data %>%
    dplyr::filter(group == "SB_negative_control") %>%
    dplyr::group_by(pos) %>%
    dplyr::summarise(coverage = mean(coverage, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(group = "SB-negative")

  sb_nonne_coverage <- coverage_data %>%
    dplyr::filter(group == "SB_positive_nonNE") %>%
    dplyr::group_by(pos) %>%
    dplyr::summarise(coverage = mean(coverage, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(group = "SB+nonNE")

  # Get sample counts for labels
  n_nepc <- length(nepc_samples)
  n_sb_neg <- length(unique(coverage_data$sample[coverage_data$group == "SB_negative_control"]))
  n_sb_nonne <- length(unique(coverage_data$sample[coverage_data$group == "SB_positive_nonNE"]))

  cat("Sample counts: NEPC=", n_nepc, ", SB-neg=", n_sb_neg, ", SB+nonNE=", n_sb_nonne, "\n")

  # Get unique insertions with orientation info
  unique_insertions <- gene_insertions %>%
    dplyr::group_by(insertion_midpoint) %>%
    dplyr::summarise(
      n_sense = sum(is_same_strand == TRUE),
      n_antisense = sum(is_same_strand == FALSE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      orientation = dplyr::case_when(
        n_sense > 0 & n_antisense > 0 ~ "Mixed",
        n_sense > 0 ~ "Sense",
        n_antisense > 0 ~ "Antisense",
        TRUE ~ "Unknown"
      )
    )

  # Get exons
  exons <- get_gene_exons(gene_symbol, gtf_df)

  # Plot range: gene boundaries + 500bp padding
  padding <- 500
  plot_range <- c(gene_start - padding, gene_end + padding)
  x_breaks <- seq(plot_range[1], plot_range[2], length.out = 5)
  x_labels <- format(round(x_breaks / 1e6, 2), nsmall = 2)

  # Filter coverage data to plot range
  nepc_coverage <- nepc_coverage %>% dplyr::filter(pos >= plot_range[1], pos <= plot_range[2])
  sb_neg_coverage <- sb_neg_coverage %>% dplyr::filter(pos >= plot_range[1], pos <= plot_range[2])
  sb_nonne_coverage <- sb_nonne_coverage %>% dplyr::filter(pos >= plot_range[1], pos <= plot_range[2])

  # Create normalized versions (scale each group to max = 1)
  nepc_coverage_norm <- nepc_coverage %>%
    dplyr::mutate(coverage = coverage / max(coverage, na.rm = TRUE))
  sb_neg_coverage_norm <- sb_neg_coverage %>%
    dplyr::mutate(coverage = coverage / max(coverage, na.rm = TRUE))
  sb_nonne_coverage_norm <- sb_nonne_coverage %>%
    dplyr::mutate(coverage = coverage / max(coverage, na.rm = TRUE))

  # Common theme
  theme_track <- theme_minimal(base_size = 11) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      plot.margin = margin(2, 10, 2, 10),
      plot.title = element_text(size = 10, face = "bold", hjust = 0)
    )

  theme_bottom <- theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_text(color = "transparent"),
      plot.margin = margin(2, 10, 5, 10)
    )

  # === Panel 4: Gene Structure with Insertions (shared by both plots) ===

  # Create gene structure data
  gene_y <- 0.25  # Y position for gene track

  # Add y-coordinates to insertion data
  unique_insertions$y_stem_bottom <- gene_y
  unique_insertions$y_stem_top <- gene_y + 0.45
  unique_insertions$y_guide_top <- 1.0

  # Color scale for orientation
  orientation_colors <- c("Sense" = COLORS$sense,
                          "Antisense" = COLORS$antisense,
                          "Mixed" = COLORS$mixed)

  # Base plot for gene structure
  p4 <- ggplot() +
    # Faint vertical guide lines extending upward
    geom_segment(data = unique_insertions,
                 aes(x = insertion_midpoint, xend = insertion_midpoint,
                     y = y_stem_top, yend = y_guide_top, color = orientation),
                 linewidth = 0.3, alpha = 0.4, linetype = "dashed") +
    # Intron line (thin line spanning gene)
    annotate("segment",
             x = gene_start, xend = gene_end,
             y = gene_y, yend = gene_y,
             color = COLORS$intron, linewidth = 0.5) +
    # Exon rectangles
    geom_rect(data = exons,
              aes(xmin = start, xmax = end, ymin = gene_y - 0.12, ymax = gene_y + 0.12),
              fill = COLORS$exon, color = COLORS$exon) +
    # Insertion lollipops - stems (colored by orientation)
    geom_segment(data = unique_insertions,
                 aes(x = insertion_midpoint, xend = insertion_midpoint,
                     y = y_stem_bottom, yend = y_stem_top, color = orientation),
                 linewidth = 1.2) +
    # Insertion lollipops - heads (colored by orientation)
    geom_point(data = unique_insertions,
               aes(x = insertion_midpoint, y = y_stem_top, fill = orientation),
               color = "white", size = 3.5, shape = 21, stroke = 0.5) +
    # Direction arrow (TSS indicator)
    annotate("segment",
             x = ifelse(gene_strand == "-", gene_end, gene_start),
             xend = ifelse(gene_strand == "-", gene_end - (gene_end - gene_start) * 0.03,
                           gene_start + (gene_end - gene_start) * 0.03),
             y = gene_y - 0.18, yend = gene_y - 0.18,
             arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
             color = COLORS$exon, linewidth = 1) +
    scale_color_manual(values = orientation_colors, name = "Orientation") +
    scale_fill_manual(values = orientation_colors, name = "Orientation") +
    scale_x_continuous(limits = plot_range, expand = c(0, 0),
                       breaks = x_breaks,
                       labels = paste0(x_labels, " Mb")) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(x = paste0("Position on ", chr), y = "Coverage",
         title = paste0(gene_symbol, " (", gene_strand, " strand)")) +
    theme_bottom +
    theme(axis.text.x = element_text(size = 9),
          plot.title = element_text(size = 10, face = "bold.italic", hjust = 0.5),
          legend.position = "inside",
          legend.position.inside = c(0.95, 0.95),
          legend.justification = c(1, 1),
          legend.background = element_rect(fill = "white", color = NA),
          legend.key.size = unit(0.4, "cm"),
          legend.title = element_text(size = 8, face = "bold"),
          legend.text = element_text(size = 7))

  # === Helper function to build combined plot ===
  build_plot <- function(nepc_cov, sb_neg_cov, sb_nonne_cov, y_label, plot_title) {
    max_cov <- max(c(nepc_cov$coverage, sb_neg_cov$coverage,
                     sb_nonne_cov$coverage), na.rm = TRUE)

    # Add vertical guide lines for insertion positions through coverage panels
    p1 <- ggplot(nepc_cov, aes(x = pos, y = coverage)) +
      geom_vline(data = unique_insertions,
                 aes(xintercept = insertion_midpoint, color = orientation),
                 linewidth = 0.4, alpha = 0.5, linetype = "dashed") +
      geom_area(fill = COLORS$nepc, alpha = 0.7, color = COLORS$nepc, linewidth = 0.3) +
      scale_color_manual(values = orientation_colors, guide = "none") +
      scale_y_continuous(limits = c(0, max_cov * 1.05), expand = c(0, 0)) +
      scale_x_continuous(limits = plot_range, expand = c(0, 0)) +
      labs(title = paste0("NEPC (n=", n_nepc, ")"), y = y_label) +
      theme_track

    p2 <- ggplot(sb_neg_cov, aes(x = pos, y = coverage)) +
      geom_vline(data = unique_insertions,
                 aes(xintercept = insertion_midpoint, color = orientation),
                 linewidth = 0.4, alpha = 0.5, linetype = "dashed") +
      geom_area(fill = COLORS$sb_neg, alpha = 0.7, color = COLORS$sb_neg, linewidth = 0.3) +
      scale_color_manual(values = orientation_colors, guide = "none") +
      scale_y_continuous(limits = c(0, max_cov * 1.05), expand = c(0, 0)) +
      scale_x_continuous(limits = plot_range, expand = c(0, 0)) +
      labs(title = paste0("SB-negative (n=", n_sb_neg, ")"), y = y_label) +
      theme_track

    p3 <- ggplot(sb_nonne_cov, aes(x = pos, y = coverage)) +
      geom_vline(data = unique_insertions,
                 aes(xintercept = insertion_midpoint, color = orientation),
                 linewidth = 0.4, alpha = 0.5, linetype = "dashed") +
      geom_area(fill = COLORS$sb_nonne, alpha = 0.7, color = COLORS$sb_nonne, linewidth = 0.3) +
      scale_color_manual(values = orientation_colors, guide = "none") +
      scale_y_continuous(limits = c(0, max_cov * 1.05), expand = c(0, 0)) +
      scale_x_continuous(limits = plot_range, expand = c(0, 0)) +
      labs(title = paste0("SB+nonNE (n=", n_sb_nonne, ")"), y = y_label) +
      theme_track

    combined <- (p1 + coord_cartesian(xlim = plot_range)) /
                (p2 + coord_cartesian(xlim = plot_range)) /
                (p3 + coord_cartesian(xlim = plot_range)) /
                (p4 + coord_cartesian(xlim = plot_range)) +
      plot_layout(heights = c(1, 1, 1, 0.8))

    combined <- combined +
      plot_annotation(
        title = plot_title,
        theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
      ) &
      theme(plot.margin = margin(2, 5, 2, 5))

    return(combined)
  }

  # === Generate RAW coverage plot ===
  raw_plot <- build_plot(
    nepc_coverage, sb_neg_coverage, sb_nonne_coverage,
    y_label = "Coverage",
    plot_title = paste0(gene_symbol, " - RNA-seq Coverage at Insertion Site")
  )

  output_file_raw <- file.path(plots_dir, paste0(gene_symbol, "_coverage.pdf"))
  ggsave(output_file_raw, raw_plot, width = 10, height = 8, device = "pdf")
  cat("Saved:", output_file_raw, "\n")

  # === Generate NORMALIZED coverage plot ===
  norm_plot <- build_plot(
    nepc_coverage_norm, sb_neg_coverage_norm, sb_nonne_coverage_norm,
    y_label = "Normalized",
    plot_title = paste0(gene_symbol, " - Normalized RNA-seq Coverage at Insertion Site")
  )

  output_file_norm <- file.path(plots_dir, paste0(gene_symbol, "_coverage_normalized.pdf"))
  ggsave(output_file_norm, norm_plot, width = 10, height = 8, device = "pdf")
  cat("Saved:", output_file_norm, "\n")

  return(c(output_file_raw, output_file_norm))
}

# === Main plotting loop ===
target_plot_genes <- c("Msi2", "Sirt1")

for (gene in target_plot_genes) {
  tryCatch({
    generate_pub_plot(gene)
  }, error = function(e) {
    cat("ERROR for", gene, ":", e$message, "\n")
  })
}

cat("\nPart 3 Complete\n")

# =============================================================================
# FINAL SUMMARY
# =============================================================================

cat("\n")
cat("=======================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("=======================================================================\n")
cat("End time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("Output Files Generated:\n")
cat("\n1. Data Preparation (Part 1):\n")
cat("   - ", output_insertions, "\n")
cat("   - ", output_samples, "\n")
cat("   - ", output_rds, "\n")

cat("\n2. Coverage Analysis (Part 2):\n")
cat("   - ", coverage_dir, "/<gene>_coverage.csv (per-gene coverage data)\n")
cat("   - ", file.path(tables_dir, "coverage_ratios.csv"), "\n")
cat("   - ", file.path(tables_dir, "analysis_summary.csv"), "\n")

cat("\n3. Publication Plots (Part 3):\n")
cat("   - ", plots_dir, "/<gene>_coverage.pdf (raw coverage)\n")
cat("   - ", plots_dir, "/<gene>_coverage_normalized.pdf (normalized)\n")

cat("\nGenes Analyzed:", paste(target_genes, collapse = ", "), "\n")
cat("Control Groups:\n")
cat("  - SB-negative controls:", length(sb_neg_controls), "samples\n")
cat("  - SB-positive nonNE controls:", length(sb_pos_nonne_controls), "samples\n")

cat("\n=======================================================================\n")
