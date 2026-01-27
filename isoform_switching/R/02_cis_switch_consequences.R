#!/usr/bin/env Rscript
#
# 02_cis_switch_consequences.R
#
# Part 2 of the CIS-integrated isoform switching analysis.
# Imports Pfam results and performs consequence analysis with protein domain information.
#
# Prerequisites:
#   - Run 01_cis_isoform_switch_overlap.R first
#   - Upload cis_genes_AA.fasta to HMMER webserver
#   - Download Pfam results file
#
# Input files:
#   - results/isoformswitch/cis_overlap/checkpoint_cis_overlap.rds
#   - Pfam results file (user-provided, from HMMER webserver)
#
# Output files:
#   - results/isoformswitch/cis_overlap/switchList_with_pfam.rds
#   - results/isoformswitch/cis_overlap/cis_switch_consequences.csv
#   - results/isoformswitch/cis_overlap/domain_switches.csv
#   - results/isoformswitch/cis_overlap/plots/*.pdf (switch plots for each gene)
#
# Usage:
#   Rscript R/02_cis_switch_consequences.R <path_to_pfam_results.txt>
#
# Example:
#   Rscript R/02_cis_switch_consequences.R results/isoformswitch/cis_overlap/pfam_results.txt
#

# === Load Required Libraries ===
suppressPackageStartupMessages({
  library(IsoformSwitchAnalyzeR)
  library(tidyverse)
})

# === Configuration ===
args <- commandArgs(trailingOnly = TRUE)

# Check for Pfam results file argument
if (length(args) < 1) {
  cat("Usage: Rscript R/02_cis_switch_consequences.R <pfam_results.txt>\n\n")
  cat("Arguments:\n")
  cat("  pfam_results.txt  Path to Pfam results file from HMMER webserver\n\n")
  cat("Example:\n")
  cat("  Rscript R/02_cis_switch_consequences.R results/isoformswitch/cis_overlap/pfam_results.txt\n\n")
  stop("Missing required argument: path to Pfam results file")
}

pfam_results_file <- args[1]

# Determine project directory
if (length(args) > 1) {
  project_dir <- args[2]
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
checkpoint_file <- file.path(project_dir, "results/isoformswitch/cis_overlap/checkpoint_cis_overlap.rds")
output_dir <- file.path(project_dir, "results/isoformswitch/cis_overlap")
plots_dir <- file.path(output_dir, "plots")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

cat("=======================================================================\n")
cat("CIS-Isoform Switch Consequence Analysis (Part 2)\n")
cat("=======================================================================\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Project directory:", project_dir, "\n")
cat("Pfam results:", pfam_results_file, "\n\n")

# === Verify Input Files ===
cat("Checking input files...\n")

if (!file.exists(checkpoint_file)) {
  stop("Checkpoint file not found: ", checkpoint_file, "\n",
       "Run 01_cis_isoform_switch_overlap.R first.")
}

if (!file.exists(pfam_results_file)) {
  stop("Pfam results file not found: ", pfam_results_file, "\n",
       "Download results from HMMER webserver first.")
}

cat("  - Checkpoint: OK\n")
cat("  - Pfam results: OK\n\n")

# === Load Checkpoint Data ===
cat("-----------------------------------------------------------------------\n")
cat("Loading Checkpoint Data\n")
cat("-----------------------------------------------------------------------\n")

checkpoint_data <- readRDS(checkpoint_file)

overlap_genes <- checkpoint_data$overlap_genes
switchList <- checkpoint_data$switchList_subset
gene_summary <- checkpoint_data$gene_summary

cat("Loaded checkpoint from:", format(checkpoint_data$analysis_date, "%Y-%m-%d %H:%M:%S"), "\n")
cat("  - Overlapping genes:", length(overlap_genes), "\n")
cat("  - Isoforms in subset:", nrow(switchList$isoformFeatures), "\n\n")

# === Import Pfam Results ===
cat("-----------------------------------------------------------------------\n")
cat("Importing Pfam Results\n")
cat("-----------------------------------------------------------------------\n")

tryCatch({
  switchList <- analyzePFAM(
    switchAnalyzeRlist = switchList,
    pathToPFAMresultFile = pfam_results_file,
    showProgress = TRUE
  )

  # Check how many domains were found
  if (!is.null(switchList$domainAnalysis) && nrow(switchList$domainAnalysis) > 0) {
    n_domains <- nrow(switchList$domainAnalysis)
    n_isoforms_with_domains <- length(unique(switchList$domainAnalysis$isoform_id))
    cat("\nPfam import successful:\n")
    cat("  - Total domain annotations:", n_domains, "\n")
    cat("  - Isoforms with domains:", n_isoforms_with_domains, "\n")

    # Show domain summary
    domain_counts <- switchList$domainAnalysis %>%
      group_by(hmm_name) %>%
      summarise(count = n(), .groups = "drop") %>%
      arrange(desc(count)) %>%
      head(10)

    cat("\nTop 10 domains found:\n")
    print(as.data.frame(domain_counts))
    cat("\n")
  } else {
    cat("\nWARNING: No protein domains found in Pfam results.\n")
    cat("This may indicate:\n")
    cat("  - Pfam results file format issue\n")
    cat("  - No domains detected in the sequences\n\n")
  }

}, error = function(e) {
  cat("\nERROR importing Pfam results:", e$message, "\n")
  cat("Please check the Pfam results file format.\n")
  cat("Expected format: Tab-delimited output from HMMER hmmscan\n\n")
  stop(e)
})

# === Analyze Alternative Splicing (prerequisite for consequences) ===
cat("-----------------------------------------------------------------------\n")
cat("Analyzing Alternative Splicing\n")
cat("-----------------------------------------------------------------------\n")

tryCatch({
  switchList <- analyzeAlternativeSplicing(
    switchList,
    quiet = FALSE
  )
  cat("Alternative splicing analysis complete.\n\n")
}, error = function(e) {
  cat("WARNING: Could not analyze alternative splicing:", e$message, "\n")
  cat("Intron retention consequences will be skipped.\n\n")
})

# === Analyze Switch Consequences (with domains) ===
cat("-----------------------------------------------------------------------\n")
cat("Analyzing Switch Consequences\n")
cat("-----------------------------------------------------------------------\n")

# Determine which consequences we can analyze based on available data
available_consequences <- c('domains_identified')

# Check if AA sequences are in switchList (needed for ORF_seq_similarity)
if (!is.null(switchList$aaSequence) && length(switchList$aaSequence) > 0) {
  available_consequences <- c(available_consequences, 'ORF_seq_similarity')
  cat("  - ORF_seq_similarity: available\n")
} else {
  cat("  - ORF_seq_similarity: skipped (no AA sequences in switchList)\n")
}

# Check if alternative splicing was analyzed (needed for intron_retention)
if (!is.null(switchList$AlternativeSplicingAnalysis)) {
  available_consequences <- c('intron_retention', available_consequences)
  cat("  - intron_retention: available\n")
} else {
  cat("  - intron_retention: skipped (no alternative splicing data)\n")
}

# Check if CPAT/CPC2 was run (needed for coding_potential)
if (!is.null(switchList$isoformFeatures$codingPotentialValue)) {
  available_consequences <- c(available_consequences, 'coding_potential')
  cat("  - coding_potential: available\n")
} else {
  cat("  - coding_potential: skipped (no CPAT/CPC2 data)\n")
}

# Check if NMD analysis was run
if (!is.null(switchList$isoformFeatures$PTC)) {
  available_consequences <- c(available_consequences, 'NMD_status')
  cat("  - NMD_status: available\n")
} else {
  cat("  - NMD_status: skipped (no NMD data)\n")
}

cat("  - domains_identified: available (from Pfam)\n\n")

tryCatch({
  switchList <- analyzeSwitchConsequences(
    switchList,
    consequencesToAnalyze = available_consequences,
    alpha = 0.05,
    dIFcutoff = 0.1,
    showProgress = TRUE
  )

  cat("\nConsequence analysis complete.\n")

  # Check consequences
  if (!is.null(switchList$switchConsequence) && nrow(switchList$switchConsequence) > 0) {
    consequence_summary <- switchList$switchConsequence %>%
      group_by(consequenceType, switchConsequence) %>%
      summarise(count = n(), .groups = "drop") %>%
      arrange(consequenceType, desc(count))

    cat("\nConsequence summary:\n")
    print(as.data.frame(consequence_summary))
    cat("\n")
  }

}, error = function(e) {
  cat("\nERROR in consequence analysis:", e$message, "\n")
  cat("Attempting minimal analysis (domains only)...\n\n")

  # Fallback to just domains
  tryCatch({
    switchList <- analyzeSwitchConsequences(
      switchList,
      consequencesToAnalyze = c('domains_identified'),
      alpha = 0.05,
      dIFcutoff = 0.1,
      showProgress = TRUE
    )
    cat("Consequence analysis (domains only) complete.\n\n")
  }, error = function(e2) {
    cat("ERROR in fallback consequence analysis:", e2$message, "\n")
  })
})

# === Extract and Save Results ===
cat("-----------------------------------------------------------------------\n")
cat("Extracting and Saving Results\n")
cat("-----------------------------------------------------------------------\n")

# 1. Save updated switchList with Pfam data
switchList_file <- file.path(output_dir, "switchList_with_pfam.rds")
saveRDS(switchList, switchList_file)
cat("Saved:", basename(switchList_file), "\n")

# 2. Extract switch consequences
if (!is.null(switchList$switchConsequence) && nrow(switchList$switchConsequence) > 0) {
  consequences_df <- switchList$switchConsequence %>%
    left_join(
      switchList$isoformFeatures %>%
        select(isoform_id, gene_name, dIF, isoform_switch_q_value),
      by = "isoform_id"
    ) %>%
    arrange(gene_name, isoform_id)

  consequences_file <- file.path(output_dir, "cis_switch_consequences.csv")
  write_csv(consequences_df, consequences_file)
  cat("Saved:", basename(consequences_file), "\n")
}

# 3. Extract domain-specific switches (if domains were found)
if (!is.null(switchList$domainAnalysis) && nrow(switchList$domainAnalysis) > 0) {
  # Get domains for switching isoforms
  switching_isoforms <- switchList$isoformFeatures %>%
    filter(isoform_switch_q_value < 0.05, abs(dIF) > 0.1) %>%
    pull(isoform_id)

  domain_switches <- switchList$domainAnalysis %>%
    filter(isoform_id %in% switching_isoforms) %>%
    left_join(
      switchList$isoformFeatures %>%
        select(isoform_id, gene_name, dIF, isoform_switch_q_value),
      by = "isoform_id"
    ) %>%
    arrange(gene_name, isoform_id)

  if (nrow(domain_switches) > 0) {
    domain_file <- file.path(output_dir, "domain_switches.csv")
    write_csv(domain_switches, domain_file)
    cat("Saved:", basename(domain_file), "\n")
  }
}

# 4. Updated gene summary with consequence info
if (!is.null(switchList$switchConsequence) && nrow(switchList$switchConsequence) > 0) {
  gene_consequences <- switchList$switchConsequence %>%
    left_join(
      switchList$isoformFeatures %>% select(isoform_id, gene_name),
      by = "isoform_id"
    ) %>%
    group_by(gene_name) %>%
    summarise(
      n_consequences = n(),
      consequence_types = paste(unique(consequenceType), collapse = "; "),
      has_domain_change = any(consequenceType == "domains_identified"),
      has_NMD_change = any(consequenceType == "NMD_status"),
      has_intron_retention = any(consequenceType == "intron_retention"),
      .groups = "drop"
    )

  # Merge with original summary
  gene_summary_updated <- gene_summary %>%
    left_join(gene_consequences, by = "gene_name")

  summary_file <- file.path(output_dir, "cis_switch_summary_with_consequences.csv")
  write_csv(gene_summary_updated, summary_file)
  cat("Saved:", basename(summary_file), "\n")
}

cat("\n")

# === Generate Switch Plots ===
cat("-----------------------------------------------------------------------\n")
cat("Generating Switch Plots for CIS Genes\n")
cat("-----------------------------------------------------------------------\n")

# Get conditions from the corrected data
# After swapping in the switchList:
#   - condition_1 = "nonNE" (reference, left bar)
#   - condition_2 = "NEPC" (comparison, right bar)
#   - IF1 = nonNE fraction, IF2 = NEPC fraction
#   - dIF = NEPC - nonNE
# So dIF > 0 means INCREASED in NEPC, dIF < 0 means DECREASED in NEPC

condition1 <- unique(switchList$isoformFeatures$condition_1)[1]  # nonNE
condition2 <- unique(switchList$isoformFeatures$condition_2)[1]  # NEPC

cat("Conditions: Reference =", condition1, ", Comparison =", condition2, "\n")
cat("Interpretation: 'Increased usage' = higher in NEPC, 'Decreased usage' = lower in NEPC\n\n")

# Generate plots for each overlapping gene
plots_generated <- 0

for (gene in overlap_genes) {
  tryCatch({
    # Create safe filename
    safe_gene_name <- gsub("[^A-Za-z0-9_-]", "_", gene)
    pdf_file <- file.path(plots_dir, paste0("switch_", safe_gene_name, ".pdf"))

    pdf(pdf_file, width = 12, height = 8, onefile = FALSE)

    switchPlot(
      switchList,
      gene = gene,
      condition1 = condition1,
      condition2 = condition2,
      localTheme = theme_bw(base_size = 12)
    )

    dev.off()
    plots_generated <- plots_generated + 1
    cat("  Generated:", basename(pdf_file), "\n")

  }, error = function(e) {
    cat("  Note: Could not create plot for gene", gene, "-", e$message, "\n")
    # Close any open device
    while (dev.cur() > 1) dev.off()
  })
}

cat("\nGenerated", plots_generated, "switch plots\n\n")

# === Print Final Summary ===
# === Generate Final Per-Gene Summary (Table S9) ===
cat("-----------------------------------------------------------------------\n")
cat("Generating Final Per-Gene Summary (Table S9)\n")
cat("-----------------------------------------------------------------------\n")

# Check if openxlsx is available for Excel export
has_openxlsx <- requireNamespace("openxlsx", quietly = TRUE)
if (has_openxlsx) {
  library(openxlsx)
}

# Load existing data files
switch_detail <- read_csv(
  file.path(output_dir, "cis_switch_genes_detail.csv"),
  show_col_types = FALSE
)

cis_insertions <- read_csv(
  file.path(output_dir, "cis_insertions_in_switching_genes.csv"),
  show_col_types = FALSE
)

# Load domain switches if available
domain_file <- file.path(output_dir, "domain_switches.csv")
has_domains <- file.exists(domain_file)

if (has_domains) {
  domain_switches <- read_csv(domain_file, show_col_types = FALSE)

  # Summarize domains per gene
  domain_summary <- domain_switches %>%
    group_by(gene_name) %>%
    summarise(
      domains_identified = paste(unique(hmm_name), collapse = "; "),
      n_unique_domains = n_distinct(hmm_name),
      n_isoforms_with_domains = n_distinct(isoform_id),
      .groups = "drop"
    )

  # Get domains per isoform with direction
  domain_by_isoform <- domain_switches %>%
    group_by(gene_name, isoform_id, dIF) %>%
    summarise(
      domains = paste(unique(hmm_name), collapse = ", "),
      .groups = "drop"
    ) %>%
    mutate(
      direction = case_when(
        dIF > 0 ~ "increased_in_NEPC",
        dIF < 0 ~ "decreased_in_NEPC",
        TRUE ~ "unchanged"
      )
    )

  # Identify domain composition changes between isoforms
  domain_changes <- domain_by_isoform %>%
    group_by(gene_name) %>%
    summarise(
      n_isoforms_in_domain_analysis = n(),
      domain_compositions = paste(
        paste0(isoform_id, " (", direction, "): ", domains),
        collapse = " | "
      ),
      unique_domain_sets = n_distinct(domains),
      .groups = "drop"
    ) %>%
    mutate(
      has_domain_composition_change = unique_domain_sets > 1
    )
}

# Summarize CIS insertions per gene
cis_summary <- cis_insertions %>%
  group_by(gene_name) %>%
  summarise(
    n_total_insertions = n(),
    n_nepc_samples = n_distinct(sample[is_nepc == TRUE]),
    n_non_nepc_samples = n_distinct(sample[is_nepc == FALSE]),
    nepc_samples = paste(unique(rnaseq_id[is_nepc == TRUE & !is.na(rnaseq_id)]), collapse = "; "),
    insertion_classifications = paste(unique(classification), collapse = "; "),
    sb_mechanisms = paste(unique(sb_mechanism), collapse = "; "),
    .groups = "drop"
  ) %>%
  mutate(
    nepc_samples = ifelse(nepc_samples == "", NA_character_, nepc_samples)
  )

# Summarize isoform switching per gene
switch_summary <- switch_detail %>%
  group_by(gene_name) %>%
  summarise(
    n_switching_isoforms = n(),
    gene_switch_q_value = min(isoform_switch_q_value),
    max_dIF = max(dIF),
    min_dIF = min(dIF),
    max_abs_dIF = max(abs(dIF)),
    top_increased_isoform = isoform_id[which.max(dIF)],
    top_increased_dIF = max(dIF),
    top_decreased_isoform = isoform_id[which.min(dIF)],
    top_decreased_dIF = min(dIF),
    all_switching_isoforms = paste(isoform_id, collapse = "; "),
    biotypes = paste(unique(iso_biotype), collapse = "; "),
    .groups = "drop"
  )

# Merge all summaries
final_summary <- switch_summary %>%
  left_join(cis_summary, by = "gene_name")

if (has_domains) {
  final_summary <- final_summary %>%
    left_join(domain_summary, by = "gene_name") %>%
    left_join(
      domain_changes %>% select(gene_name, has_domain_composition_change, domain_compositions),
      by = "gene_name"
    ) %>%
    mutate(
      has_domain_annotation = !is.na(domains_identified),
      has_domain_composition_change = replace_na(has_domain_composition_change, FALSE)
    )
} else {
  final_summary <- final_summary %>%
    mutate(
      has_domain_annotation = FALSE,
      n_unique_domains = NA_integer_,
      domains_identified = NA_character_,
      has_domain_composition_change = FALSE,
      domain_compositions = NA_character_
    )
}

# Reorder and select columns
final_summary <- final_summary %>%
  select(
    gene_name,
    n_switching_isoforms,
    gene_switch_q_value,
    max_abs_dIF,
    top_increased_isoform,
    top_increased_dIF,
    top_decreased_isoform,
    top_decreased_dIF,
    biotypes,
    n_total_insertions,
    n_nepc_samples,
    n_non_nepc_samples,
    nepc_samples,
    insertion_classifications,
    sb_mechanisms,
    has_domain_annotation,
    n_unique_domains,
    domains_identified,
    has_domain_composition_change,
    domain_compositions,
    all_switching_isoforms
  ) %>%
  arrange(gene_switch_q_value)

# Create clean manuscript table (Table S9)
table_s9 <- final_summary %>%
  select(
    Gene = gene_name,
    `N Switching Isoforms` = n_switching_isoforms,
    `Q-value` = gene_switch_q_value,
    `Max |dIF|` = max_abs_dIF,
    `Top Increased Isoform` = top_increased_isoform,
    `dIF (Increased)` = top_increased_dIF,
    `Top Decreased Isoform` = top_decreased_isoform,
    `dIF (Decreased)` = top_decreased_dIF,
    `Biotypes` = biotypes,
    `N CIS Insertions` = n_total_insertions,
    `N NEPC Samples` = n_nepc_samples,
    `NEPC Samples` = nepc_samples,
    `Insertion Types` = insertion_classifications,
    `SB Mechanism` = sb_mechanisms,
    `Has Domain Annotation` = has_domain_annotation,
    `N Domains` = n_unique_domains,
    `Domains` = domains_identified
  )

# Save CSV files
write_csv(final_summary, file.path(output_dir, "cis_isoform_switch_summary_full.csv"))
cat("Saved: cis_isoform_switch_summary_full.csv\n")

write_csv(table_s9, file.path(output_dir, "Table_S9_isoform_switching.csv"))
cat("Saved: Table_S9_isoform_switching.csv\n")

# Save Excel file with formatting
if (has_openxlsx) {
  wb <- createWorkbook()

  # Sheet 1: Main summary table (Table S9)
  addWorksheet(wb, "Table S9")
  writeData(wb, "Table S9", table_s9)

  # Format header
  headerStyle <- createStyle(
    textDecoration = "bold",
    fgFill = "#4472C4",
    fontColour = "white",
    halign = "center",
    border = "TopBottomLeftRight"
  )
  addStyle(wb, "Table S9", headerStyle, rows = 1, cols = 1:ncol(table_s9))

  # Format data cells
  dataStyle <- createStyle(border = "TopBottomLeftRight", halign = "left")
  addStyle(wb, "Table S9", dataStyle,
           rows = 2:(nrow(table_s9)+1),
           cols = 1:ncol(table_s9),
           gridExpand = TRUE)

  # Number formatting for q-values
  qvalStyle <- createStyle(numFmt = "0.00E+00", halign = "right")
  addStyle(wb, "Table S9", qvalStyle,
           rows = 2:(nrow(table_s9)+1),
           cols = 3,
           gridExpand = TRUE)

  # Number formatting for dIF
  difStyle <- createStyle(numFmt = "0.00", halign = "right")
  addStyle(wb, "Table S9", difStyle,
           rows = 2:(nrow(table_s9)+1),
           cols = c(4, 6, 8),
           gridExpand = TRUE)

  setColWidths(wb, "Table S9", cols = 1:ncol(table_s9), widths = "auto")

  # Sheet 2: Domain details (if available)
  if (has_domains) {
    addWorksheet(wb, "Domain Details")
    writeData(wb, "Domain Details", domain_by_isoform)
    addStyle(wb, "Domain Details", headerStyle, rows = 1, cols = 1:ncol(domain_by_isoform))
    setColWidths(wb, "Domain Details", cols = 1:ncol(domain_by_isoform), widths = "auto")
  }

  # Sheet 3: Full data
  addWorksheet(wb, "Full Data")
  writeData(wb, "Full Data", final_summary)
  addStyle(wb, "Full Data", headerStyle, rows = 1, cols = 1:ncol(final_summary))
  setColWidths(wb, "Full Data", cols = 1:ncol(final_summary), widths = "auto")

  saveWorkbook(wb, file.path(output_dir, "Table_S9_isoform_switching.xlsx"), overwrite = TRUE)
  cat("Saved: Table_S9_isoform_switching.xlsx\n")
} else {
  cat("Note: openxlsx not available, skipping Excel export\n")
}

# Print summary statistics
cat("\nFinal Summary Statistics:\n")
cat("  - Total genes:", nrow(final_summary), "\n")
cat("  - Genes with domain annotation:", sum(final_summary$has_domain_annotation), "\n")
cat("  - Genes with domain composition changes:", sum(final_summary$has_domain_composition_change), "\n")

cat("\nTop 5 genes by |dIF|:\n")
final_summary %>%
  arrange(desc(max_abs_dIF)) %>%
  head(5) %>%
  select(gene_name, max_abs_dIF, n_nepc_samples, domains_identified) %>%
  print()

cat("\n")

cat("=======================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("=======================================================================\n")
cat("End time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("Summary:\n")
cat("  - Genes analyzed:", length(overlap_genes), "\n")
cat("  - Switch plots generated:", plots_generated, "\n")

if (!is.null(switchList$domainAnalysis) && nrow(switchList$domainAnalysis) > 0) {
  cat("  - Protein domains found:", nrow(switchList$domainAnalysis), "\n")
}

if (!is.null(switchList$switchConsequence) && nrow(switchList$switchConsequence) > 0) {
  cat("  - Switch consequences:", nrow(switchList$switchConsequence), "\n")
}

cat("\nOutput files:\n")
cat("  Directory:", output_dir, "\n")

output_files <- list.files(output_dir, pattern = "\\.(csv|rds)$")
for (f in output_files) {
  cat("  -", f, "\n")
}

cat("\nPlot files:\n")
cat("  Directory:", plots_dir, "\n")

plot_files <- list.files(plots_dir, pattern = "\\.pdf$")
for (f in head(plot_files, 10)) {
  cat("  -", f, "\n")
}
if (length(plot_files) > 10) {
  cat("  ... and", length(plot_files) - 10, "more\n")
}

cat("\n=======================================================================\n")
cat("CIS-integrated isoform switching analysis complete!\n")
cat("=======================================================================\n")
