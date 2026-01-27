#!/usr/bin/env Rscript
#
# 01_cis_isoform_switch_overlap.R
#
# Integrates CIS (Common Insertion Site) data with isoform switching analysis.
# Identifies genes with both transposon insertions AND significant isoform switches.
#
# This is Part 1 of the CIS-integrated analysis workflow:
#   Part 1: Identify overlapping genes, extract sequences for Pfam analysis
#   Part 2: Import Pfam results and complete consequence analysis
#
# Input files:
#   - results/isoformswitch/checkpoints/checkpoint_NEPC_vs_nonNE_03_dexseq.rds
#   - ../fusion_finder/R/all_insertions_master.csv (all CIS genes)
#   - ../fusion_finder/R/key_genes_insertions.csv (key CIS genes)
#
# Output files:
#   - results/isoformswitch/cis_overlap/cis_switch_overlap_summary.csv
#   - results/isoformswitch/cis_overlap/cis_switch_genes_detail.csv
#   - results/isoformswitch/cis_overlap/cis_switch_sample_patterns.csv
#   - results/isoformswitch/cis_overlap/sequences/cis_genes_AA.fasta (for HMMER upload)
#   - results/isoformswitch/cis_overlap/checkpoint_cis_overlap.rds
#
# Usage: Rscript R/01_cis_isoform_switch_overlap.R
#
# After running this script:
#   1. Upload cis_genes_AA.fasta to HMMER webserver (https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan)
#   2. Download Pfam results
#   3. Run R/02_cis_switch_consequences.R with the Pfam results
#

# === Load Required Libraries ===
suppressPackageStartupMessages({
  library(IsoformSwitchAnalyzeR)
  library(tidyverse)
  library(Biostrings)  # For sequence handling
})

# === Configuration ===
args <- commandArgs(trailingOnly = TRUE)

if (length(args) > 0) {
  project_dir <- args[1]
} else {
  if (file.exists("data/sample_metadata.csv")) {
    project_dir <- "."
  } else if (file.exists("../data/sample_metadata.csv")) {
    project_dir <- ".."
  } else {
    stop("Cannot find project directory. Run from isoform_switching/ directory.")
  }
}

# Define paths
checkpoint_file <- file.path(project_dir, "results/isoformswitch/checkpoints/checkpoint_NEPC_vs_nonNE_03_dexseq.rds")
all_insertions_file <- file.path(project_dir, "../fusion_finder/R/all_insertions_master.csv")
key_genes_file <- file.path(project_dir, "../fusion_finder/R/key_genes_insertions.csv")

# Output directory
output_dir <- file.path(project_dir, "results/isoformswitch/cis_overlap")
sequences_dir <- file.path(output_dir, "sequences")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(sequences_dir, recursive = TRUE, showWarnings = FALSE)

cat("=======================================================================\n")
cat("CIS-Isoform Switch Overlap Analysis (Part 1)\n")
cat("=======================================================================\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Project directory:", project_dir, "\n\n")

# === Verify Input Files ===
cat("Checking input files...\n")

if (!file.exists(checkpoint_file)) {
  stop("Checkpoint file not found: ", checkpoint_file, "\n",
       "Run isoform_switch_analysis.R first to generate DEXSeq results.")
}

if (!file.exists(all_insertions_file)) {
  stop("CIS insertions file not found: ", all_insertions_file)
}

cat("  - Checkpoint: OK\n")
cat("  - CIS insertions: OK\n\n")

# === Load Data ===
cat("-----------------------------------------------------------------------\n")
cat("Loading Data\n")
cat("-----------------------------------------------------------------------\n")

# Load switchAnalyzeRlist
cat("Loading isoform switch data...\n")
switchList <- readRDS(checkpoint_file)

cat("  - Isoforms:", nrow(switchList$isoformFeatures), "\n")
cat("  - Genes:", length(unique(switchList$isoformFeatures$gene_name)), "\n")
cat("  - Samples:", nrow(switchList$designMatrix), "\n\n")

# Load CIS insertion data
cat("Loading CIS insertion data...\n")
all_insertions <- read_csv(all_insertions_file, show_col_types = FALSE)
cat("  - Total insertions:", nrow(all_insertions), "\n")
cat("  - Unique CIS genes:", length(unique(all_insertions$gene_name)), "\n")

# Load key genes if available
if (file.exists(key_genes_file)) {
  key_genes <- read_csv(key_genes_file, show_col_types = FALSE)
  key_gene_names <- unique(key_genes$gene_name)
  cat("  - Key CIS genes:", length(key_gene_names), "\n")
} else {
  key_gene_names <- NULL
  cat("  - Key genes file not found (optional)\n")
}

cat("\n")

# === Extract Significant Isoform Switches ===
cat("-----------------------------------------------------------------------\n")
cat("Identifying Significant Isoform Switches\n")
cat("-----------------------------------------------------------------------\n")

# Get significant switches (q < 0.05, |dIF| > 0.1)
sig_switches <- switchList$isoformFeatures %>%
  filter(
    isoform_switch_q_value < 0.05,
    abs(dIF) > 0.1
  )

switching_genes <- unique(sig_switches$gene_name)

cat("Significant switches:\n")
cat("  - Isoforms:", nrow(sig_switches), "\n")
cat("  - Genes:", length(switching_genes), "\n\n")

# === Find Overlapping Genes ===
cat("-----------------------------------------------------------------------\n")
cat("Finding CIS-Switch Overlapping Genes\n")
cat("-----------------------------------------------------------------------\n")

# Get unique CIS genes
cis_genes <- unique(all_insertions$gene_name)

# Find overlap
overlap_genes <- intersect(switching_genes, cis_genes)

cat("Overlap analysis:\n")
cat("  - CIS genes:", length(cis_genes), "\n")
cat("  - Switching genes:", length(switching_genes), "\n")
cat("  - Overlapping genes:", length(overlap_genes), "\n")

if (length(overlap_genes) == 0) {
  cat("\nWARNING: No overlapping genes found!\n")
  cat("This may indicate:\n")
  cat("  - Gene name mismatches between datasets\n")
  cat("  - CIS genes don't show isoform switching in this comparison\n\n")
} else {
  cat("\nOverlapping genes:\n")
  cat(paste("  ", overlap_genes, collapse = "\n"), "\n\n")
}

# Check key genes overlap
if (!is.null(key_gene_names)) {
  key_overlap <- intersect(switching_genes, key_gene_names)
  cat("Key gene overlap:\n")
  cat("  - Key CIS genes with switches:", length(key_overlap), "\n")
  if (length(key_overlap) > 0) {
    cat("  - Genes:", paste(key_overlap, collapse = ", "), "\n")
  }
  cat("\n")
}

# === Detailed Analysis of Overlapping Genes ===
cat("-----------------------------------------------------------------------\n")
cat("Analyzing Overlapping Genes in Detail\n")
cat("-----------------------------------------------------------------------\n")

if (length(overlap_genes) > 0) {

  # Get switch details for overlapping genes
  overlap_switches <- sig_switches %>%
    filter(gene_name %in% overlap_genes) %>%
    arrange(isoform_switch_q_value)

  cat("Switch details for overlapping genes:\n")
  cat("  - Isoforms involved:", nrow(overlap_switches), "\n")

  # Get CIS insertion details for overlapping genes
  overlap_insertions <- all_insertions %>%
    filter(gene_name %in% overlap_genes) %>%
    select(gene_name, sample, rnaseq_id, is_nepc, classification,
           sb_mechanism, insertion_chr, insertion_start, insertion_end) %>%
    distinct()

  cat("  - CIS insertions in these genes:", nrow(overlap_insertions), "\n\n")

  # === Sample-Level Analysis ===
  # For each overlapping gene, check if samples with insertions show different patterns

  cat("Analyzing sample-level isoform expression patterns...\n\n")

  # Get sample-level expression data
  iso_expr <- switchList$isoformRepExpression
  design <- switchList$designMatrix

  sample_patterns <- list()

  for (gene in overlap_genes) {
    # Get isoforms for this gene
    gene_isoforms <- switchList$isoformFeatures %>%
      filter(gene_name == gene) %>%
      pull(isoform_id)

    # Get samples with insertions in this gene
    gene_insertions <- overlap_insertions %>%
      filter(gene_name == gene, !is.na(rnaseq_id))

    # Parse comma-separated rnaseq_ids
    insertion_samples <- gene_insertions %>%
      pull(rnaseq_id) %>%
      str_split(",") %>%
      unlist() %>%
      str_trim() %>%
      unique() %>%
      na.omit()

    # Get expression for this gene's isoforms across samples
    if (length(gene_isoforms) > 0 && any(gene_isoforms %in% iso_expr$isoform_id)) {
      gene_expr <- iso_expr %>%
        filter(isoform_id %in% gene_isoforms) %>%
        pivot_longer(cols = -isoform_id, names_to = "sample", values_to = "expression") %>%
        left_join(design, by = c("sample" = "sampleID")) %>%
        mutate(
          has_insertion = sample %in% insertion_samples,
          gene_name = gene
        )

      sample_patterns[[gene]] <- gene_expr
    }
  }

  # Combine all sample patterns
  if (length(sample_patterns) > 0) {
    all_sample_patterns <- bind_rows(sample_patterns)

    # Summarize by gene
    pattern_summary <- all_sample_patterns %>%
      group_by(gene_name, isoform_id, condition, has_insertion) %>%
      summarise(
        mean_expr = mean(expression, na.rm = TRUE),
        sd_expr = sd(expression, na.rm = TRUE),
        n_samples = n(),
        .groups = "drop"
      )

    cat("Sample pattern analysis complete.\n")
    cat("  - Genes analyzed:", length(unique(pattern_summary$gene_name)), "\n\n")
  }

  # === Create Summary Tables ===
  cat("-----------------------------------------------------------------------\n")
  cat("Creating Summary Tables\n")
  cat("-----------------------------------------------------------------------\n")

  # 1. Gene-level overlap summary
  gene_summary <- overlap_switches %>%
    group_by(gene_name) %>%
    summarise(
      n_switching_isoforms = n(),
      min_q_value = min(isoform_switch_q_value),
      max_abs_dIF = max(abs(dIF)),
      isoforms = paste(isoform_id, collapse = "; "),
      .groups = "drop"
    ) %>%
    left_join(
      overlap_insertions %>%
        group_by(gene_name) %>%
        summarise(
          n_insertions = n(),
          n_samples_with_insertion = n_distinct(rnaseq_id, na.rm = TRUE),
          insertion_classifications = paste(unique(classification), collapse = "; "),
          sb_mechanisms = paste(unique(sb_mechanism), collapse = "; "),
          .groups = "drop"
        ),
      by = "gene_name"
    ) %>%
    mutate(
      is_key_gene = gene_name %in% key_gene_names
    ) %>%
    arrange(min_q_value)

  # Save gene summary
  summary_file <- file.path(output_dir, "cis_switch_overlap_summary.csv")
  write_csv(gene_summary, summary_file)
  cat("Saved:", basename(summary_file), "\n")

  # 2. Detailed isoform-level results
  detail_df <- overlap_switches %>%
    select(gene_name, isoform_id, IF1, IF2, dIF,
           isoform_switch_q_value, gene_switch_q_value,
           iso_biotype) %>%
    left_join(
      overlap_insertions %>%
        group_by(gene_name) %>%
        summarise(
          cis_samples = paste(unique(na.omit(rnaseq_id)), collapse = "; "),
          cis_classifications = paste(unique(classification), collapse = "; "),
          .groups = "drop"
        ),
      by = "gene_name"
    ) %>%
    arrange(isoform_switch_q_value)

  detail_file <- file.path(output_dir, "cis_switch_genes_detail.csv")
  write_csv(detail_df, detail_file)
  cat("Saved:", basename(detail_file), "\n")

  # 3. Sample-level patterns
  if (exists("all_sample_patterns") && nrow(all_sample_patterns) > 0) {
    patterns_file <- file.path(output_dir, "cis_switch_sample_patterns.csv")
    write_csv(all_sample_patterns, patterns_file)
    cat("Saved:", basename(patterns_file), "\n")
  }

  # 4. CIS insertion details for overlapping genes
  insertions_file <- file.path(output_dir, "cis_insertions_in_switching_genes.csv")
  write_csv(overlap_insertions, insertions_file)
  cat("Saved:", basename(insertions_file), "\n")

  cat("\n")

  # === Extract Amino Acid Sequences for Pfam Analysis ===
  cat("-----------------------------------------------------------------------\n")
  cat("Extracting Amino Acid Sequences for Pfam Analysis\n")
  cat("-----------------------------------------------------------------------\n")

  # Get isoform IDs for overlapping genes
  overlap_isoform_ids <- switchList$isoformFeatures %>%
    filter(gene_name %in% overlap_genes) %>%
    pull(isoform_id) %>%
    unique()

  cat("Isoforms to extract:", length(overlap_isoform_ids), "\n")

  # Extract sequences using IsoformSwitchAnalyzeR
  # First, subset the switchList to only overlapping genes for efficiency
  switchList_subset <- switchList
  switchList_subset$isoformFeatures <- switchList_subset$isoformFeatures %>%
    filter(gene_name %in% overlap_genes)

  # Extract sequences
  tryCatch({
    cat("Running extractSequence()...\n")

    switchList_subset <- extractSequence(
      switchList_subset,
      pathToOutput = sequences_dir,
      writeSequences = TRUE,
      outputPrefix = "cis_genes",
      extractNTseq = FALSE,  # Only need AA for Pfam
      extractAAseq = TRUE,
      removeShortAAseq = TRUE,
      removeLongAAseq = FALSE,
      quiet = FALSE
    )

    # Check output
    aa_file <- file.path(sequences_dir, "cis_genes_AA.fasta")
    if (file.exists(aa_file)) {
      # Count sequences
      aa_seqs <- readAAStringSet(aa_file)
      cat("\nAmino acid sequences extracted:\n")
      cat("  - File:", basename(aa_file), "\n")
      cat("  - Sequences:", length(aa_seqs), "\n")
      cat("  - Ready for HMMER/Pfam upload\n")
    } else {
      cat("\nWARNING: AA FASTA file was not created.\n")
      cat("This may happen if no ORFs were found in the overlapping genes.\n")
    }

  }, error = function(e) {
    cat("\nERROR extracting sequences:", e$message, "\n")
    cat("You may need to run extractSequence manually.\n")
  })

  cat("\n")

  # === Save Checkpoint ===
  cat("-----------------------------------------------------------------------\n")
  cat("Saving Checkpoint\n")
  cat("-----------------------------------------------------------------------\n")

  checkpoint_data <- list(
    overlap_genes = overlap_genes,
    key_overlap_genes = if (!is.null(key_gene_names)) key_overlap else NULL,
    overlap_switches = overlap_switches,
    overlap_insertions = overlap_insertions,
    gene_summary = gene_summary,
    switchList_subset = switchList_subset,
    analysis_date = Sys.time()
  )

  checkpoint_output <- file.path(output_dir, "checkpoint_cis_overlap.rds")
  saveRDS(checkpoint_data, checkpoint_output)
  cat("Saved checkpoint:", basename(checkpoint_output), "\n\n")

} else {
  cat("No overlapping genes found. Skipping detailed analysis.\n\n")
}

# === Print Summary ===
cat("=======================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("=======================================================================\n")
cat("End time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("Summary:\n")
cat("  - CIS genes analyzed:", length(cis_genes), "\n")
cat("  - Switching genes:", length(switching_genes), "\n")
cat("  - Overlapping genes:", length(overlap_genes), "\n")
if (!is.null(key_gene_names)) {
  cat("  - Key genes with switches:", length(key_overlap), "\n")
}

cat("\nOutput files:\n")
output_files <- list.files(output_dir, recursive = TRUE)
for (f in output_files) {
  cat("  -", f, "\n")
}

cat("\n")
cat("=======================================================================\n")
cat("NEXT STEPS\n")
cat("=======================================================================\n")
cat("1. Upload the AA FASTA file to HMMER webserver:\n")
cat("   https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan\n\n")
cat("2. Select 'Pfam' as the target database\n\n")
cat("3. Download the results (tab-delimited format)\n\n")
cat("4. Run Part 2 of the analysis:\n")
cat("   Rscript R/02_cis_switch_consequences.R <pfam_results.txt>\n")
cat("=======================================================================\n")
