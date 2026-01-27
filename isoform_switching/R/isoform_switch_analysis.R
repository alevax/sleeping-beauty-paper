#!/usr/bin/env Rscript
#
# isoform_switch_analysis.R
#
# Comprehensive Isoform Switching Analysis using IsoformSwitchAnalyzeR
#
# Comparisons:
#   1. NEPC vs nonNE (within SBT2+)
#   2. SBT2+ (NEPC + nonNE) vs SBT2_neg
#
# This script performs:
#   - Data import from Salmon quantifications
#   - Isoform switch detection with DEXSeq
#   - ORF prediction and analysis
#   - Functional consequence prediction
#   - Comprehensive visualizations
#   - Summary statistics and reports
#
# Prerequisites:
#   - Salmon quantification complete (results/salmon/)
#   - Sample metadata (data/sample_metadata.csv)
#   - GENCODE GTF annotation (reference/gencode.vM25.annotation.gtf)
#   - GENCODE transcriptome FASTA (reference/gencode.vM25.transcripts.fa)
#
# Usage: Rscript R/isoform_switch_analysis.R
#

# === Load Required Libraries ===
suppressPackageStartupMessages({
  library(IsoformSwitchAnalyzeR)
  library(tidyverse)
  library(tximport)  # For importing Salmon quantifications
  library(BSgenome)  # For sequence analysis if needed
})

# === Fix dplyr compatibility issue ===
# IsoformSwitchAnalyzeR uses old dplyr syntax with 'all=TRUE' in full_join
# This has been deprecated in favor of keeping all rows by default
# We need to create a wrapper that ignores the 'all' parameter
if (packageVersion("dplyr") >= "1.0.0") {
  cat("Applying dplyr compatibility fix for IsoformSwitchAnalyzeR...\n")

  # Store original full_join
  full_join_original <- dplyr::full_join

  # Create wrapper that removes 'all' parameter
  full_join_fixed <- function(x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"), ..., all = NULL) {
    # Ignore 'all' parameter and call original without it
    full_join_original(x, y, by = by, copy = copy, suffix = suffix, ...)
  }

  # Replace in dplyr namespace
  assignInNamespace("full_join", full_join_fixed, ns = "dplyr")

  cat("dplyr compatibility fix applied successfully.\n\n")
}

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
    stop("Cannot find sample_metadata.csv. Run from isoform_switching/ directory.")
  }
}

# Define paths
salmon_dir <- file.path(project_dir, "results/salmon")
metadata_file <- file.path(project_dir, "data/sample_metadata.csv")
gtf_file <- file.path(project_dir, "reference/gencode.vM25.annotation.gtf")
fasta_file <- file.path(project_dir, "reference/gencode.vM25.transcripts.fa")
output_dir <- file.path(project_dir, "results/isoformswitch")
plots_dir <- file.path(output_dir, "plots")

# Create output directories
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

# Create checkpoint directory for saving intermediate results
checkpoint_dir <- file.path(output_dir, "checkpoints")
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

# === Checkpoint Helper Functions ===
# These functions enable automatic recovery from long-running analyses
# Checkpoints are saved after expensive operations (import, DEXSeq, consequences)
# If a checkpoint exists, the step is skipped and data is loaded from disk

save_checkpoint <- function(obj, comparison_name, step_name) {
  checkpoint_file <- file.path(checkpoint_dir,
    paste0("checkpoint_", comparison_name, "_", step_name, ".rds"))
  saveRDS(obj, checkpoint_file)
  cat("  [Checkpoint saved:", basename(checkpoint_file), "]\n")
}

load_checkpoint <- function(comparison_name, step_name) {
  checkpoint_file <- file.path(checkpoint_dir,
    paste0("checkpoint_", comparison_name, "_", step_name, ".rds"))
  if (file.exists(checkpoint_file)) {
    cat("  [Loading checkpoint:", basename(checkpoint_file), "]\n")
    return(readRDS(checkpoint_file))
  }
  return(NULL)
}

checkpoint_exists <- function(comparison_name, step_name) {
  checkpoint_file <- file.path(checkpoint_dir,
    paste0("checkpoint_", comparison_name, "_", step_name, ".rds"))
  file.exists(checkpoint_file)
}

# Setup logging
log_file <- file.path(output_dir, paste0("analysis_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message")

cat("=======================================================================\n")
cat("IsoformSwitchAnalyzeR - Comprehensive Isoform Switching Analysis\n")
cat("=======================================================================\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Project directory:", project_dir, "\n")
cat("Salmon dir:", salmon_dir, "\n")
cat("Metadata:", metadata_file, "\n")
cat("GTF file:", gtf_file, "\n")
cat("FASTA file:", fasta_file, "\n")
cat("Output:", output_dir, "\n\n")

# Verify files exist
if (!file.exists(metadata_file)) stop("Metadata file not found!")
if (!file.exists(gtf_file)) stop("GTF file not found!")
if (!file.exists(fasta_file)) stop("FASTA file not found!")
if (!dir.exists(salmon_dir)) stop("Salmon directory not found!")

# === Load Sample Metadata ===
cat("-----------------------------------------------------------------------\n")
cat("Loading Sample Metadata\n")
cat("-----------------------------------------------------------------------\n")

metadata <- read_csv(metadata_file, show_col_types = FALSE)

# Filter to known samples only (exclude 'unknown' condition if any)
metadata <- metadata %>%
  filter(condition %in% c("NEPC", "nonNE", "SBT2_neg"))

cat("Total samples loaded:", nrow(metadata), "\n")
cat("\nSamples by condition:\n")
print(table(metadata$condition))
cat("\n")

# Check which samples have Salmon output
metadata <- metadata %>%
  mutate(
    salmon_path = file.path(salmon_dir, sample),
    has_salmon = dir.exists(salmon_path) &
                 file.exists(file.path(salmon_path, "quant.sf"))
  )

cat("Samples with complete Salmon output:", sum(metadata$has_salmon), "/", nrow(metadata), "\n")

if (sum(metadata$has_salmon) == 0) {
  stop("No Salmon output found. Run Salmon quantification first.")
}

# Filter to samples with Salmon output
metadata <- metadata %>% filter(has_salmon)

cat("\nFinal samples by condition:\n")
condition_table <- table(metadata$condition)
print(condition_table)
cat("\n")

# === Verify Salmon Output Files ===
cat("-----------------------------------------------------------------------\n")
cat("Verifying Salmon Output Files\n")
cat("-----------------------------------------------------------------------\n")

# Verify all quant.sf files exist
all_samples <- metadata$sample
for (sample_id in all_samples) {
  quant_file <- file.path(salmon_dir, sample_id, "quant.sf")
  if (!file.exists(quant_file)) {
    stop("Missing quant.sf for sample: ", sample_id)
  }
}

cat("All", length(all_samples), "Salmon quant.sf files verified successfully.\n")
cat("Ready to run comparisons.\n\n")

# =============================================================================
# HELPER FUNCTIONS FOR COMPREHENSIVE ANALYSIS
# =============================================================================

#' Run comprehensive switching analysis for a comparison
#'
#' @param comparison_name Name of the comparison (for output files)
#' @param design_matrix Data frame with sampleID and condition columns
#' @param condition1 First condition name
#' @param condition2 Second condition name
#' @return switchAnalyzeRlist object with complete analysis
run_comprehensive_analysis <- function(comparison_name, design_matrix, condition1, condition2) {

  cat("=======================================================================\n")
  cat("COMPARISON:", comparison_name, "\n")
  cat("=======================================================================\n")
  cat("Conditions:", condition1, "vs", condition2, "\n")
  cat("Sample sizes:",
      sum(design_matrix$condition == condition1), "(", condition1, ") vs",
      sum(design_matrix$condition == condition2), "(", condition2, ")\n\n")

  # --- Step 1: Import Salmon data and Create switchAnalyzeRlist ---
  # CHECKPOINT: After import (saves 5-15 minutes on re-run)
  cat("Step 1: Importing Salmon quantifications and creating switchAnalyzeRlist...\n")

  if (checkpoint_exists(comparison_name, "01_imported")) {
    switchList <- load_checkpoint(comparison_name, "01_imported")
    cat("  Skipping import - loaded from checkpoint\n")
    cat("  switchList contains:\n")
    cat("    - Isoforms:", nrow(switchList$isoformFeatures), "\n")
    cat("    - Genes:", length(unique(switchList$isoformFeatures$gene_id)), "\n\n")
  } else {
    tryCatch({
      # Use tximport to import Salmon quantification files
      salmon_files <- file.path(salmon_dir, design_matrix$sampleID, "quant.sf")
      names(salmon_files) <- design_matrix$sampleID

      cat("  Importing quantification data with tximport...\n")
      txi <- tximport(salmon_files, type = "salmon", txOut = TRUE)

      cat("  Successfully imported data for", ncol(txi$abundance), "samples and",
          nrow(txi$abundance), "transcripts\n")

      # Convert to data frames with isoform_id column for IsoformSwitchAnalyzeR
      abundance_df <- as.data.frame(txi$abundance)
      abundance_df$isoform_id <- rownames(abundance_df)
      abundance_df <- abundance_df[, c("isoform_id", design_matrix$sampleID)]

      counts_df <- as.data.frame(txi$counts)
      counts_df$isoform_id <- rownames(counts_df)
      counts_df <- counts_df[, c("isoform_id", design_matrix$sampleID)]

      cat("  Creating switchAnalyzeRlist with importRdata...\n")

      # Due to dplyr compatibility issues, we need to temporarily use older join behavior
      # Create the switchList without GTF first, then add it separately
      switchList <- suppressWarnings(
        importRdata(
          isoformCountMatrix = counts_df,
          isoformRepExpression = abundance_df,
          designMatrix = design_matrix,
          isoformExonAnnoation = gtf_file,
          isoformNtFasta = fasta_file,
          showProgress = TRUE,
          quiet = FALSE,
          addAnnotatedORFs = TRUE,
          removeTECgenes = TRUE
        )
      )

      cat("  Created switchList with:\n")
      cat("    - Isoforms:", nrow(switchList$isoformFeatures), "\n")
      cat("    - Genes:", length(unique(switchList$isoformFeatures$gene_id)), "\n\n")

      # Save checkpoint after successful import
      save_checkpoint(switchList, comparison_name, "01_imported")

    }, error = function(e) {
      cat("ERROR in data import:", e$message, "\n")
      cat("Full error trace:\n")
      print(e)
      stop(e)
    })
  }

  # --- Step 2: Filter low expression ---
  cat("Step 2: Filtering lowly expressed isoforms...\n")

  pre_filter_isoforms <- nrow(switchList$isoformFeatures)
  pre_filter_genes <- length(unique(switchList$isoformFeatures$gene_id))

  switchList <- preFilter(
    switchList,
    geneExpressionCutoff = 1,
    isoformExpressionCutoff = 0,
    IFcutoff = 0.01,  # Minimum isoform fraction
    removeSingleIsoformGenes = TRUE,
    quiet = FALSE
  )

  post_filter_isoforms <- nrow(switchList$isoformFeatures)
  post_filter_genes <- length(unique(switchList$isoformFeatures$gene_id))

  cat("  Filtering results:\n")
  cat("    - Isoforms: ", pre_filter_isoforms, "->", post_filter_isoforms,
      "(", round(100*post_filter_isoforms/pre_filter_isoforms, 1), "%)\n")
  cat("    - Genes: ", pre_filter_genes, "->", post_filter_genes,
      "(", round(100*post_filter_genes/pre_filter_genes, 1), "%)\n\n")

  # --- Step 3: Test for differential isoform usage ---
  # CHECKPOINT: After DEXSeq (saves 10-30 minutes on re-run)
  cat("Step 3: Testing for differential isoform usage with DEXSeq...\n")

  if (checkpoint_exists(comparison_name, "03_dexseq")) {
    switchList <- load_checkpoint(comparison_name, "03_dexseq")
    cat("  Skipping DEXSeq - loaded from checkpoint\n\n")
  } else {
    switchList <- isoformSwitchTestDEXSeq(
      switchAnalyzeRlist = switchList,
      reduceToSwitchingGenes = FALSE,  # Keep all for now
      reduceFurtherToGenesWithConsequencePotential = FALSE,
      onlySigIsoforms = FALSE,
      keepIsoformInAllConditions = TRUE,
      alpha = 0.05,
      dIFcutoff = 0.1  # Minimum isoform fraction change
    )

    cat("  DEXSeq testing complete.\n")

    # Save checkpoint after DEXSeq
    save_checkpoint(switchList, comparison_name, "03_dexseq")
    cat("\n")
  }

  # --- Step 4: Extract ORF information ---
  cat("Step 4: Analyzing Open Reading Frames (ORFs)...\n")

  tryCatch({
    switchList <- analyzeORF(
      switchList,
      genomeObject = NULL,  # Use FASTA sequences already loaded
      orfMethod = "longest",
      minORFlength = 30,  # Minimum 10 amino acids
      showProgress = TRUE
    )

    cat("  ORF analysis complete.\n")
    cat("    - Isoforms with ORF predicted:",
        sum(!is.na(switchList$orfAnalysis$orfTransciptStart)), "\n\n")

  }, error = function(e) {
    cat("  WARNING: ORF analysis failed:", e$message, "\n")
    cat("  Continuing without ORF information...\n\n")
  })

  # --- Step 5: Analyze Alternative Splicing ---
  # IMPORTANT: This must be run BEFORE analyzeSwitchConsequences() to enable
  # intron_retention analysis. This step classifies alternative splicing events.
  cat("Step 5: Analyzing alternative splicing events...\n")

  tryCatch({
    switchList <- analyzeAlternativeSplicing(
      switchList,
      onlySwitchingGenes = TRUE,
      alpha = 0.05,
      dIFcutoff = 0.1,
      showProgress = TRUE,
      quiet = FALSE
    )

    cat("  Alternative splicing analysis complete.\n")

    # Report splicing event counts
    if ("AlternativeSplicingAnalysis" %in% names(switchList) &&
        !is.null(switchList$AlternativeSplicingAnalysis)) {
      as_analysis <- switchList$AlternativeSplicingAnalysis

      cat("  Splicing events detected:\n")
      if ("IR" %in% colnames(as_analysis)) {
        cat("    - Intron Retention (IR):", sum(as_analysis$IR > 0, na.rm = TRUE), "isoforms\n")
      }
      if ("ES" %in% colnames(as_analysis)) {
        cat("    - Exon Skipping (ES):", sum(as_analysis$ES > 0, na.rm = TRUE), "isoforms\n")
      }
      if ("MEE" %in% colnames(as_analysis)) {
        cat("    - Mutually Exclusive Exons (MEE):", sum(as_analysis$MEE > 0, na.rm = TRUE), "isoforms\n")
      }
      if ("A5" %in% colnames(as_analysis)) {
        cat("    - Alternative 5' Splice Site (A5):", sum(as_analysis$A5 > 0, na.rm = TRUE), "isoforms\n")
      }
      if ("A3" %in% colnames(as_analysis)) {
        cat("    - Alternative 3' Splice Site (A3):", sum(as_analysis$A3 > 0, na.rm = TRUE), "isoforms\n")
      }
      if ("ATSS" %in% colnames(as_analysis)) {
        cat("    - Alternative TSS (ATSS):", sum(as_analysis$ATSS > 0, na.rm = TRUE), "isoforms\n")
      }
      if ("ATTS" %in% colnames(as_analysis)) {
        cat("    - Alternative TTS (ATTS):", sum(as_analysis$ATTS > 0, na.rm = TRUE), "isoforms\n")
      }
    }
    cat("\n")

  }, error = function(e) {
    cat("  WARNING: Alternative splicing analysis failed:", e$message, "\n")
    cat("  Continuing without splicing classification...\n\n")
  })

  # --- Step 6: Analyze switching consequences ---
  # CHECKPOINT: After consequences (saves 20-50 minutes on re-run)
  cat("Step 6: Analyzing functional switching consequences...\n")

  if (checkpoint_exists(comparison_name, "06_consequences")) {
    switchList <- load_checkpoint(comparison_name, "06_consequences")
    cat("  Skipping consequence analysis - loaded from checkpoint\n\n")
  } else {
    tryCatch({
      switchList <- analyzeSwitchConsequences(
        switchList,
        consequencesToAnalyze = c(
          'intron_retention',
          'ORF_seq_similarity',
          'NMD_status'
        ),
        alpha = 0.05,
        dIFcutoff = 0.1,
        showProgress = TRUE
      )

      cat("  Consequence analysis complete.\n")

      # Save checkpoint after all analysis steps complete
      save_checkpoint(switchList, comparison_name, "06_consequences")
      cat("\n")

    }, error = function(e) {
      cat("  WARNING: Consequence analysis failed:", e$message, "\n")
      cat("  Continuing without consequence information...\n\n")
    })
  }

  # --- Step 7: Extract and filter to significant switches ---
  cat("Step 7: Extracting significant isoform switches...\n")

  # Extract summary before filtering
  full_summary <- extractSwitchSummary(switchList, alpha = 0.05, dIFcutoff = 0.1)

  cat("  Summary of all genes:\n")
  print(full_summary)
  cat("\n")

  # Filter to switching genes
  switchList <- extractTopSwitches(
    switchList,
    filterForConsequences = FALSE,
    n = Inf,  # Keep all significant
    sortByQvals = TRUE
  )

  # Get switching summary
  switch_summary <- extractSwitchSummary(switchList, alpha = 0.05, dIFcutoff = 0.1)

  cat("  Significant isoform switches detected:\n")
  cat("    - Switching genes:", switch_summary$nrGenes[1], "\n")
  cat("    - Switching isoforms:", switch_summary$nrIsoforms[1], "\n\n")

  # --- Step 8: Save results ---
  cat("Step 8: Saving results...\n")

  # Save RDS object
  rds_file <- file.path(output_dir, paste0("switchList_", comparison_name, ".rds"))
  saveRDS(switchList, rds_file)
  cat("  Saved: switchList_", comparison_name, ".rds\n", sep = "")

  # Export switching results if any found
  if (switch_summary$nrGenes[1] > 0) {

    # Extract top switches
    top_switches <- extractTopSwitches(
      switchList,
      filterForConsequences = FALSE,
      n = Inf,
      sortByQvals = TRUE,
      extractGenes = FALSE
    )

    switches_file <- file.path(output_dir, paste0("top_switches_", comparison_name, ".csv"))
    write_csv(as.data.frame(top_switches), switches_file)
    cat("  Saved: top_switches_", comparison_name, ".csv (",
        nrow(top_switches), " isoforms)\n", sep = "")

    # Extract gene-level summary
    gene_summary <- extractTopSwitches(
      switchList,
      filterForConsequences = FALSE,
      n = Inf,
      sortByQvals = TRUE,
      extractGenes = TRUE
    )

    genes_file <- file.path(output_dir, paste0("switching_genes_", comparison_name, ".csv"))
    write_csv(as.data.frame(gene_summary), genes_file)
    cat("  Saved: switching_genes_", comparison_name, ".csv (",
        nrow(gene_summary), " genes)\n", sep = "")

    # Extract consequence enrichment if available
    if ("switchConsequence" %in% names(switchList)) {
      tryCatch({
        consequence_enrich <- extractConsequenceEnrichment(
          switchList,
          alpha = 0.05,
          dIFcutoff = 0.1,
          countGenes = TRUE
        )

        if (!is.null(consequence_enrich) && nrow(consequence_enrich) > 0) {
          enrich_file <- file.path(output_dir, paste0("consequence_enrichment_", comparison_name, ".csv"))
          write_csv(consequence_enrich, enrich_file)
          cat("  Saved: consequence_enrichment_", comparison_name, ".csv\n", sep = "")
        }
      }, error = function(e) {
        cat("  Note: Could not extract consequence enrichment\n")
      })
    }

    # Extract splicing summary
    tryCatch({
      splicing_summary <- extractSplicingSummary(
        switchList,
        alpha = 0.05,
        dIFcutoff = 0.1,
        onlySigIsoforms = FALSE
      )

      if (!is.null(splicing_summary) && nrow(splicing_summary) > 0) {
        splicing_file <- file.path(output_dir, paste0("splicing_summary_", comparison_name, ".csv"))
        write_csv(splicing_summary, splicing_file)
        cat("  Saved: splicing_summary_", comparison_name, ".csv\n", sep = "")
      }
    }, error = function(e) {
      cat("  Note: Could not extract splicing summary\n")
    })

  } else {
    cat("  No significant switches found - skipping detailed exports\n")
  }

  cat("\n")

  # --- Step 9: Generate visualizations ---
  if (switch_summary$nrGenes[1] > 0) {
    cat("Step 9: Generating visualizations...\n")

    # Switch consequence summary plot
    tryCatch({
      pdf_file <- file.path(plots_dir, paste0("consequence_summary_", comparison_name, ".pdf"))
      pdf(pdf_file, width = 10, height = 8)

      switchPlotTopSwitches(
        switchList,
        filterForConsequences = FALSE,
        n = min(20, switch_summary$nrGenes[1]),
        sortByQvals = TRUE
      )

      dev.off()
      cat("  Saved: consequence_summary_", comparison_name, ".pdf\n", sep = "")
    }, error = function(e) {
      cat("  Note: Could not create top switches plot:", e$message, "\n")
    })

    # Volcano-like plot
    tryCatch({
      pdf_file <- file.path(plots_dir, paste0("volcano_plot_", comparison_name, ".pdf"))
      pdf(pdf_file, width = 10, height = 8)

      switchPlot(
        switchList,
        gene = NULL,  # Overview plot
        condition1 = condition1,
        condition2 = condition2
      )

      dev.off()
      cat("  Saved: volcano_plot_", comparison_name, ".pdf\n", sep = "")
    }, error = function(e) {
      cat("  Note: Could not create volcano plot:", e$message, "\n")
    })

    # Individual switch plots for top genes (max 10)
    top_genes <- unique(top_switches$gene_name)[1:min(10, switch_summary$nrGenes[1])]

    for (gene in top_genes) {
      tryCatch({
        safe_gene_name <- gsub("[^A-Za-z0-9_-]", "_", gene)
        pdf_file <- file.path(plots_dir, paste0("switch_", comparison_name, "_", safe_gene_name, ".pdf"))

        pdf(pdf_file, width = 12, height = 8)

        switchPlot(
          switchList,
          gene = gene,
          condition1 = condition1,
          condition2 = condition2,
          localTheme = theme_bw()
        )

        dev.off()
      }, error = function(e) {
        cat("  Note: Could not create plot for gene", gene, "\n")
      })
    }

    cat("  Saved individual gene plots for top", length(top_genes), "genes\n")
    cat("\n")
  }

  # --- Step 10: Global Alternative Splicing Analysis ---
  cat("Step 10: Performing global alternative splicing analysis...\n")

  tryCatch({
    # Extract splicing summary
    splicing_summary_global <- extractSplicingSummary(
      switchList,
      alpha = 0.05,
      dIFcutoff = 0.1,
      onlySigIsoforms = FALSE,
      plot = FALSE,
      returnResult = TRUE
    )

    if (!is.null(splicing_summary_global) && nrow(splicing_summary_global) > 0) {
      splicing_global_file <- file.path(output_dir, paste0("splicing_summary_global_", comparison_name, ".csv"))
      write_csv(splicing_summary_global, splicing_global_file)
      cat("  Saved: splicing_summary_global_", comparison_name, ".csv\n", sep = "")
    }

    # Splicing enrichment analysis
    tryCatch({
      pdf_file <- file.path(plots_dir, paste0("splicing_enrichment_", comparison_name, ".pdf"))
      pdf(pdf_file, width = 10, height = 8)

      splicing_enrichment <- extractSplicingEnrichment(
        switchList,
        alpha = 0.05,
        dIFcutoff = 0.1,
        onlySigIsoforms = FALSE,
        returnResult = TRUE,
        returnSummary = FALSE
      )

      dev.off()
      cat("  Saved: splicing_enrichment_", comparison_name, ".pdf\n", sep = "")

      if (!is.null(splicing_enrichment) && nrow(splicing_enrichment) > 0) {
        splicing_enrich_file <- file.path(output_dir, paste0("splicing_enrichment_", comparison_name, ".csv"))
        write_csv(splicing_enrichment, splicing_enrich_file)
        cat("  Saved: splicing_enrichment_", comparison_name, ".csv\n", sep = "")
      }
    }, error = function(e) {
      cat("  Note: Could not generate splicing enrichment:", e$message, "\n")
    })

    # Genome-wide splicing analysis
    tryCatch({
      pdf_file <- file.path(plots_dir, paste0("splicing_genomewide_", comparison_name, ".pdf"))
      pdf(pdf_file, width = 12, height = 8)

      splicing_genomewide <- extractSplicingGenomeWide(
        switchList,
        alpha = 0.05,
        dIFcutoff = 0.1,
        onlySigIsoforms = FALSE,
        plot = TRUE,
        returnResult = TRUE
      )

      dev.off()
      cat("  Saved: splicing_genomewide_", comparison_name, ".pdf\n", sep = "")

      if (!is.null(splicing_genomewide) && nrow(splicing_genomewide) > 0) {
        splicing_gw_file <- file.path(output_dir, paste0("splicing_genomewide_", comparison_name, ".csv"))
        write_csv(splicing_genomewide, splicing_gw_file)
        cat("  Saved: splicing_genomewide_", comparison_name, ".csv\n", sep = "")
      }
    }, error = function(e) {
      cat("  Note: Could not generate genome-wide splicing analysis:", e$message, "\n")
    })

    cat("\n")

  }, error = function(e) {
    cat("  WARNING: Global splicing analysis failed:", e$message, "\n\n")
  })

  # --- Step 11: Global Consequence Analysis ---
  cat("Step 11: Performing global consequence analysis...\n")

  tryCatch({
    # Consequence summary
    tryCatch({
      pdf_file <- file.path(plots_dir, paste0("consequence_summary_global_", comparison_name, ".pdf"))
      pdf(pdf_file, width = 12, height = 8)

      consequence_summary <- extractConsequenceSummary(
        switchList,
        consequencesToAnalyze = 'all',
        alpha = 0.05,
        dIFcutoff = 0.1,
        plot = TRUE,
        plotGenes = FALSE,
        asFractionTotal = FALSE,
        returnResult = TRUE
      )

      dev.off()
      cat("  Saved: consequence_summary_global_", comparison_name, ".pdf\n", sep = "")

      if (!is.null(consequence_summary) && nrow(consequence_summary) > 0) {
        conseq_summary_file <- file.path(output_dir, paste0("consequence_summary_global_", comparison_name, ".csv"))
        write_csv(consequence_summary, conseq_summary_file)
        cat("  Saved: consequence_summary_global_", comparison_name, ".csv\n", sep = "")
      }
    }, error = function(e) {
      cat("  Note: Could not generate consequence summary:", e$message, "\n")
    })

    # Consequence enrichment analysis (gain vs loss)
    tryCatch({
      pdf_file <- file.path(plots_dir, paste0("consequence_enrichment_global_", comparison_name, ".pdf"))
      pdf(pdf_file, width = 12, height = 8)

      consequence_enrichment <- extractConsequenceEnrichment(
        switchList,
        consequencesToAnalyze = 'all',
        alpha = 0.05,
        dIFcutoff = 0.1,
        analysisOppositeConsequence = TRUE,
        plot = TRUE,
        returnResult = TRUE
      )

      dev.off()
      cat("  Saved: consequence_enrichment_global_", comparison_name, ".pdf\n", sep = "")

      if (!is.null(consequence_enrichment) && nrow(consequence_enrichment) > 0) {
        conseq_enrich_file <- file.path(output_dir, paste0("consequence_enrichment_global_", comparison_name, ".csv"))
        write_csv(consequence_enrichment, conseq_enrich_file)
        cat("  Saved: consequence_enrichment_global_", comparison_name, ".csv\n", sep = "")
      }
    }, error = function(e) {
      cat("  Note: Could not generate consequence enrichment:", e$message, "\n")
    })

    # Genome-wide consequence analysis
    tryCatch({
      pdf_file <- file.path(plots_dir, paste0("consequence_genomewide_", comparison_name, ".pdf"))
      pdf(pdf_file, width = 12, height = 10)

      consequence_genomewide <- extractConsequenceGenomeWide(
        switchList,
        consequencesToAnalyze = 'all',
        alpha = 0.05,
        dIFcutoff = 0.1,
        plot = TRUE,
        returnResult = TRUE
      )

      dev.off()
      cat("  Saved: consequence_genomewide_", comparison_name, ".pdf\n", sep = "")

      if (!is.null(consequence_genomewide) && nrow(consequence_genomewide) > 0) {
        conseq_gw_file <- file.path(output_dir, paste0("consequence_genomewide_", comparison_name, ".csv"))
        write_csv(consequence_genomewide, conseq_gw_file)
        cat("  Saved: consequence_genomewide_", comparison_name, ".csv\n", sep = "")
      }
    }, error = function(e) {
      cat("  Note: Could not generate genome-wide consequence analysis:", e$message, "\n")
    })

    cat("\n")

  }, error = function(e) {
    cat("  WARNING: Global consequence analysis failed:", e$message, "\n\n")
  })

  cat("-----------------------------------------------------------------------\n")
  cat("Comparison", comparison_name, "complete!\n")
  cat("-----------------------------------------------------------------------\n\n")

  return(switchList)
}

# =============================================================================
# COMPARISON 1: NEPC vs nonNE (within SBT2+)
# =============================================================================

nepc_vs_nonne_samples <- metadata %>%
  filter(condition %in% c("NEPC", "nonNE"))

if (nrow(nepc_vs_nonne_samples) >= 6 &&
    sum(nepc_vs_nonne_samples$condition == "NEPC") >= 3 &&
    sum(nepc_vs_nonne_samples$condition == "nonNE") >= 3) {

  design_nepc <- nepc_vs_nonne_samples %>%
    select(sampleID = sample, condition) %>%
    as.data.frame()

  switchList_nepc <- run_comprehensive_analysis(
    comparison_name = "NEPC_vs_nonNE",
    design_matrix = design_nepc,
    condition1 = "nonNE",
    condition2 = "NEPC"
  )

} else {
  cat("=======================================================================\n")
  cat("SKIPPING: NEPC vs nonNE comparison\n")
  cat("=======================================================================\n")
  cat("Insufficient samples for robust analysis.\n")
  cat("Need at least 3 samples per group.\n\n")
}

# =============================================================================
# COMPARISON 2: SBT2+ vs SBT2_neg
# =============================================================================

sbt2_samples <- metadata %>%
  mutate(sbt2_group = ifelse(condition == "SBT2_neg", "SBT2_neg", "SBT2_pos"))

if (nrow(sbt2_samples) >= 6 &&
    sum(sbt2_samples$sbt2_group == "SBT2_neg") >= 3 &&
    sum(sbt2_samples$sbt2_group == "SBT2_pos") >= 3) {

  design_sbt2 <- sbt2_samples %>%
    select(sampleID = sample, condition = sbt2_group) %>%
    as.data.frame()

  switchList_sbt2 <- run_comprehensive_analysis(
    comparison_name = "SBT2pos_vs_SBT2neg",
    design_matrix = design_sbt2,
    condition1 = "SBT2_pos",
    condition2 = "SBT2_neg"
  )

} else {
  cat("=======================================================================\n")
  cat("SKIPPING: SBT2+ vs SBT2- comparison\n")
  cat("=======================================================================\n")
  cat("Insufficient samples for robust analysis.\n")
  cat("Need at least 3 samples per group.\n\n")
}

# =============================================================================
# FINAL SUMMARY
# =============================================================================

cat("\n")
cat("=======================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("=======================================================================\n")
cat("End time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Results directory:", output_dir, "\n")
cat("Plots directory:", plots_dir, "\n")
cat("Log file:", log_file, "\n")
cat("\n")

# List all output files
output_files <- list.files(output_dir, full.names = FALSE, recursive = FALSE)
cat("Output files generated:\n")
for (f in output_files) {
  cat("  -", f, "\n")
}

cat("\nPlot files generated:\n")
plot_files <- list.files(plots_dir, full.names = FALSE)
for (f in plot_files) {
  cat("  -", f, "\n")
}

cat("\n")
cat("Please review the results and visualizations.\n")
cat("For detailed switch information, see the CSV files in:", output_dir, "\n")
cat("For visualizations, see the PDF files in:", plots_dir, "\n")
cat("=======================================================================\n")

# Close log file
sink(type = "message")
sink(type = "output")
close(log_con)

cat("\nLog saved to:", log_file, "\n")
