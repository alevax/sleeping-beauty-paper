# Updated Lollipop Plot Script
# Improvements: exonic/intronic annotation, cleaner aesthetics, proper legends

# Output directories
base.dir <- "experiments/lollipop"
plots.dir <- file.path(base.dir, "plots")
tables.dir <- file.path(base.dir, "tables")

# Key genes for filtered output (reviewer-requested)
KEY_GENES <- c("Ezh2", "Sirt1", "Arid1a", "Aifm3", "Rad21", "Msi2", "Cxxc5",
               "Golph3", "Gdi2", "Prkci", "Nisch", "Spred1", "Tgfbr2",
               "Ppp2r2a", "Evi2a")

suppressPackageStartupMessages({
  library(trackViewer)
  library(GenomicRanges)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(dplyr)
  library(RColorBrewer)
})

# Load data once
load_lollipop_data <- function() {
  cat("Loading data...\n")
  CIS_.001 <<- read.csv("experiments/table/TAPDANCE_001.csv")
  txdb <<- TxDb.Mmusculus.UCSC.mm10.knownGene
  all.exonic_parts <<- suppressMessages(exonicParts(txdb, linked.to.single.gene.only=TRUE))
  all.gene_id <<- mcols(all.exonic_parts)$gene_id
  all.genes <<- suppressMessages(genes(txdb))
  all.promoters <<- promoters(all.genes, upstream=100, downstream=50)
  simplified_bed <<- read.table("TAPDANCE_NEPC/results/all/all-nr-prostate-0.001.BED",
                                 sep = "\t", header = FALSE, skip = 1)

  # Load sample metadata for NEPC status and RNA-seq mapping
  cat("Loading sample metadata...\n")
  sample_metadata <<- load_sample_metadata()
  cat("Data loaded.\n")
}

#' Load sample metadata with NEPC status and RNA-seq mapping
#' Creates a lookup table from TAPDANCE sample ID to phenotype and RNA-seq status
load_sample_metadata <- function() {
  suppressPackageStartupMessages(library(readxl))

  # Load TAPDANCE metadata (DNA sample IDs)
  tapdance <- read.csv("data/tapdance/mouse_metadata_20220912_FA.csv",
                       stringsAsFactors = FALSE)
  tapdance$mouse_id <- gsub(" .*", "", tapdance$Tumor.ID)

  # Load RNA-seq metadata from Excel - use sheet 5 "RNA-Seq Mod by Ale" which has QC-passed samples
  rnaseq <- tryCatch({
    sheet5 <- read_excel("data/mouse-analysis/SB List CIS and RNAseq_2020_OCT_20.xlsx",
                         sheet = "RNA-Seq Mod by Ale")
    colnames(sheet5) <- make.names(colnames(sheet5))
    sheet5 %>%
      filter(!is.na(RNAseq.ID) & grepl("CMZ", RNAseq.ID)) %>%
      select(RNAseq.ID, mouse.ID)
  }, error = function(e) {
    cat("Warning: Could not load RNA-seq metadata\n")
    data.frame(RNAseq.ID = character(), mouse.ID = character())
  })

  # Aggregate RNA-seq IDs for mice with multiple samples (e.g., mouse 3267 has CMZ111 and CMZ356)
  rnaseq_agg <- rnaseq %>%
    group_by(mouse.ID) %>%
    summarize(RNAseq.ID = paste(sort(RNAseq.ID), collapse = ","), .groups = "drop")

  # Create mapping: TAPDANCE sample -> has_rnaseq
  tapdance$has_rnaseq <- tapdance$mouse_id %in% rnaseq_agg$mouse.ID

  # Add RNAseq ID where available (now aggregated, so no duplicate rows)
  tapdance <- merge(tapdance, rnaseq_agg,
                    by.x = "mouse_id", by.y = "mouse.ID",
                    all.x = TRUE)

  # Create final lookup table
  metadata <- data.frame(
    sample_id = tapdance$DNA.Sample.ID,
    mouse_id = tapdance$mouse_id,
    tumor_id = tapdance$Tumor.ID,
    is_nepc = tapdance$NE.Phenotype,
    is_castrated = tapdance$Castration,
    is_metastasis = tapdance$Metastasis.sample,
    has_rnaseq = tapdance$has_rnaseq,
    rnaseq_id = tapdance$RNAseq.ID,
    stringsAsFactors = FALSE
  )

  cat("  Loaded metadata for", nrow(metadata), "samples\n")
  cat("  NEPC samples:", sum(metadata$is_nepc), "\n")
  cat("  Samples with RNA-seq:", sum(metadata$has_rnaseq), "\n")

  return(metadata)
}

# Color palette - professional and colorblind-friendly
COLORS <- list(
  exonic_same = "#27ae60",      # Green - exonic, same strand
  exonic_opposite = "#e74c3c",  # Red - exonic, opposite strand
  intronic_same = "#3498db",    # Blue - intronic, same strand
  intronic_opposite = "#9b59b6", # Purple - intronic, opposite strand
  exon_fill = "#bdc3c7",        # Light gray for exons
  exon_border = "#7f8c8d",      # Darker gray for exon borders
  promoter_fill = "#f39c12",    # Orange for promoter
  promoter_border = "#d68910"   # Darker orange for promoter border
)

# Updated lollipop function
lolliplot_updated <- function(gene_symbol, gene_entrez, promoter_bp = 5000,
                               overhang_bp = 0, show_labels = FALSE) {

  cat("\n=== Generating lollipop for:", gene_symbol, "(Entrez:", gene_entrez, ") ===\n")

  # Get exons for this gene
  gene_exons <- all.exonic_parts[all.gene_id %in% gene_entrez]
  cat("Exons found:", length(gene_exons), "\n")

  if(length(gene_exons) == 0) {
    cat("ERROR: No exons found for this gene\n")
    return(NULL)
  }

  # Style exons (NULL names to suppress legend entirely)
  gene_exons$fill <- COLORS$exon_fill
  gene_exons$color <- COLORS$exon_border
  gene_exons$height <- 0.08  # Thicker exons
  names(gene_exons) <- NULL  # NULL = no legend entries

  # Get strand and chromosome
  strand <- names(sort(table(strand(gene_exons)), decreasing=TRUE)[1])
  chromosome <- names(sort(table(seqnames(gene_exons)), decreasing=TRUE)[1])

  # Get TSS position for marker
  promoter <- all.promoters[all.promoters$gene_id %in% gene_entrez]
  cat("Promoters found:", length(promoter), "\n")

  # Create TSS marker as a small GRanges (will be drawn as vertical bar)
  # TSS is at the 5' end of the gene - where transcription starts
  if(length(gene_exons) > 0) {
    if(strand == "+") {
      tss_pos <- min(start(gene_exons))  # For + strand, TSS at lowest coord
    } else {
      tss_pos <- max(end(gene_exons))    # For - strand, TSS at highest coord
    }
    tss_marker <- GRanges(seqnames = chromosome,
                          ranges = IRanges(start = tss_pos, width = 50),  # Slightly wider TSS marker
                          strand = strand)
    tss_marker$fill <- COLORS$promoter_fill
    tss_marker$color <- COLORS$promoter_border
    tss_marker$height <- 0.12  # Thicker than exons (0.08)
    names(tss_marker) <- NULL
  } else {
    tss_marker <- NULL
  }
  cat("Strand:", strand, "| Chromosome:", chromosome, "\n")

  # Calculate gene boundaries
  if(strand == "+") {
    gene_start <- min(start(gene_exons)) - promoter_bp
    gene_end <- max(end(gene_exons)) + overhang_bp
  } else {
    gene_start <- max(end(gene_exons)) + promoter_bp
    gene_end <- min(start(gene_exons)) - overhang_bp
  }
  start_bp <- min(gene_start, gene_end)
  end_bp <- max(gene_start, gene_end)

  # Create range for plotting (tighter bounds - no extra padding)
  range <- GRanges(seqnames = chromosome,
                   ranges = IRanges(start = start_bp,
                                    end = end_bp),
                   strand = strand)

  # Get CIS data for this gene
  CIS <- CIS_.001 %>% filter(gene_name == gene_symbol)
  cat("CIS entries:", nrow(CIS), "\n")

  if(nrow(CIS) == 0) {
    cat("ERROR: No CIS entries found\n")
    return(NULL)
  }

  CIS_samples <- trimws(unlist(strsplit(CIS$library_name, "::")))
  CIS_samples <- CIS_samples[CIS_samples != ""]
  cat("Samples with insertions:", length(CIS_samples), "\n")

  # Subset insertions to gene region
  CIS_subset <- simplified_bed %>%
    filter(V1 == chromosome) %>%
    filter(V3 > start_bp) %>%
    filter(V2 < end_bp) %>%
    mutate(midpoint = (V2 + V3) / 2) %>%
    filter(V4 %in% CIS_samples)

  cat("Insertions in region:", nrow(CIS_subset), "\n")

  if(nrow(CIS_subset) == 0) {
    cat("WARNING: No insertions found in region\n")
    return(NULL)
  }

  # Create insertions GRanges
  insertions.gr <- GRanges(
    seqnames = rep(chromosome, nrow(CIS_subset)),
    ranges = IRanges(CIS_subset$midpoint, width = 1)
  )

  # Determine if each insertion overlaps exons (exonic vs intronic)
  is_exonic <- countOverlaps(insertions.gr, gene_exons) > 0

  # Determine strand relationship
  insertion_strand <- CIS_subset$V6
  if(strand == "+") {
    is_same_strand <- insertion_strand == "+"
  } else {
    is_same_strand <- insertion_strand == "-"
  }

  # Assign colors based on exonic status AND strand
  insertions.gr$color <- case_when(
    is_exonic & is_same_strand ~ COLORS$exonic_same,
    is_exonic & !is_same_strand ~ COLORS$exonic_opposite,
    !is_exonic & is_same_strand ~ COLORS$intronic_same,
    !is_exonic & !is_same_strand ~ COLORS$intronic_opposite
  )

  # Store classification for export
  insertions.gr$is_exonic <- is_exonic
  insertions.gr$is_same_strand <- is_same_strand
  insertions.gr$classification <- case_when(
    is_exonic & is_same_strand ~ "exonic_same",
    is_exonic & !is_same_strand ~ "exonic_opposite",
    !is_exonic & is_same_strand ~ "intronic_same",
    !is_exonic & !is_same_strand ~ "intronic_opposite"
  )

  # Style insertions
  insertions.gr$score <- CIS_subset$V5
  insertions.gr$height <- CIS_subset$V5
  insertions.gr$border <- "#2c3e50"
  insertions.gr$alpha <- 0.85
  insertions.gr$cex <- 0.8  # Larger dots
  insertions.gr$shape <- "circle"
  insertions.gr$SNPsideID <- "top"

  if(show_labels) {
    insertions.gr$label <- CIS_subset$V4
    insertions.gr$label.parameter.gp <- gpar(cex = 0.4, col = "#2c3e50")
    insertions.gr$label.parameter.rot <- 45
  }

  # Print summary
  cat("\n--- Insertion Classification ---\n")
  cat("Exonic (same strand):", sum(is_exonic & is_same_strand), "\n")
  cat("Exonic (opposite strand):", sum(is_exonic & !is_same_strand), "\n")
  cat("Intronic (same strand):", sum(!is_exonic & is_same_strand), "\n")
  cat("Intronic (opposite strand):", sum(!is_exonic & !is_same_strand), "\n")

  # X-axis labels
  xaxis <- floor(seq(start_bp, end_bp, length.out = 6))
  names(xaxis) <- format(xaxis, big.mark = ",", scientific = FALSE)

  # Create legend
  legend_list <- list(
    labels = c("Exonic (same strand)", "Exonic (opposite strand)",
               "Intronic (same strand)", "Intronic (opposite strand)"),
    col = c(COLORS$exonic_same, COLORS$exonic_opposite,
            COLORS$intronic_same, COLORS$intronic_opposite),
    fill = c(COLORS$exonic_same, COLORS$exonic_opposite,
             COLORS$intronic_same, COLORS$intronic_opposite)
  )

  # Ensure output directory exists
  dir.create(plots.dir, showWarnings = FALSE, recursive = TRUE)

  # Generate plot (compact dimensions)
  output_file <- file.path(plots.dir, paste0(gene_symbol, ".pdf"))
  pdf(output_file, width = 10, height = 4)

  # Combine exons with TSS marker
  if(!is.null(tss_marker)) {
    features <- c(gene_exons, tss_marker)
  } else {
    features <- gene_exons
  }

  lolliplot(insertions.gr,
            features,
            range,
            ylab = "Read Count",
            ylab.gp = gpar(fontsize = 11, fontface = "bold"),
            yaxis.gp = gpar(fontsize = 9),
            xaxis.gp = gpar(fontsize = 8),
            xaxis = xaxis,
            legend = FALSE,  # Disable auto-legend entirely, we draw our own
            jitter = "node",  # node = spread overlapping dots
            type = "circle")

  # Add custom legend (moved down further)
  pushViewport(viewport(x = 0.5, y = 0.84, width = 0.9, height = 0.05))
  grid.draw(legendGrob(
    labels = c("Exonic (same strand)", "Exonic (opposite strand)",
               "Intronic (same strand)", "Intronic (opposite strand)"),
    pch = 21,
    gp = gpar(col = c(COLORS$exonic_same, COLORS$exonic_opposite,
                      COLORS$intronic_same, COLORS$intronic_opposite),
              fill = c(COLORS$exonic_same, COLORS$exonic_opposite,
                       COLORS$intronic_same, COLORS$intronic_opposite)),
    ncol = 4,
    byrow = TRUE
  ))
  popViewport()

  # Add title (larger)
  grid.text(paste0(gene_symbol, " (", strand, " strand)"),
            x = 0.5, y = 0.91,
            gp = gpar(cex = 1.4, fontface = "bold.italic"))

  dev.off()
  cat("\nPlot saved to:", output_file, "\n")

  # Export insertion data with classifications
  export_df <- data.frame(
    sample = CIS_subset$V4,
    chromosome = CIS_subset$V1,
    start = CIS_subset$V2,
    end = CIS_subset$V3,
    midpoint = CIS_subset$midpoint,
    read_count = CIS_subset$V5,
    insertion_strand = CIS_subset$V6,
    gene_strand = strand,
    is_exonic = is_exonic,
    is_same_strand = is_same_strand,
    classification = insertions.gr$classification
  )

  # No longer saving per-gene CSVs - all data goes to master table
  return(export_df)
}

# ============================================================================
# BATCH PROCESSING FUNCTIONS
# ============================================================================

#' Run lollipop analysis for all CIS genes and build comprehensive table
#' @param distance_threshold Max distance from gene body for exonic/intronic classification (default 2000bp)
#' @param generate_plots Whether to generate PDF plots (set FALSE for faster table-only run)
#' @return Comprehensive data frame with one row per insertion
run_full_analysis <- function(distance_threshold = 2000, generate_plots = TRUE) {

  cat("========================================\n")
  cat("COMPREHENSIVE LOLLIPOP ANALYSIS\n")
  cat("========================================\n")
  cat("Distance threshold:", distance_threshold, "bp\n")
  cat("Insertions >", distance_threshold, "bp from gene body classified as 'distant'\n\n")

  # Ensure data is loaded
  if(!exists("CIS_.001")) {
    load_lollipop_data()
  }

  # Get unique gene names
  all_genes <- unique(CIS_.001$gene_name)
  all_genes <- all_genes[all_genes != "No results within 20000 bp window"]
  cat("Found", length(all_genes), "CIS genes to process\n")

  # Use biomaRt for gene symbol to Entrez ID conversion
  cat("Converting gene symbols to Entrez IDs using biomaRt...\n")
  suppressPackageStartupMessages(library(biomaRt))
  ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

  conversion <- getBM(
    attributes = c('mgi_symbol', 'entrezgene_id'),
    filters = 'mgi_symbol',
    values = all_genes,
    mart = ensembl
  )

  # Remove NA Entrez IDs and duplicates (keep first)
  conversion <- conversion[!is.na(conversion$entrezgene_id), ]
  conversion <- conversion[!duplicated(conversion$mgi_symbol), ]
  cat("Successfully mapped", nrow(conversion), "genes to Entrez IDs\n")

  # Track results
  all_insertions <- list()
  failed_genes <- c()
  skipped_genes <- setdiff(all_genes, conversion$mgi_symbol)

  cat("\nSkipped genes (no Entrez ID):", length(skipped_genes), "\n")

  # Process each gene
  cat("\n=== Processing", nrow(conversion), "genes ===\n")

  for(i in 1:nrow(conversion)) {
    gene_symbol <- conversion$mgi_symbol[i]
    gene_entrez <- as.character(conversion$entrezgene_id[i])

    cat("[", i, "/", nrow(conversion), "]", gene_symbol, "...")

    result <- tryCatch({
      process_gene(gene_symbol, gene_entrez, distance_threshold, generate_plots)
    }, error = function(e) {
      cat(" ERROR:", e$message, "\n")
      return(NULL)
    })

    if(!is.null(result) && nrow(result) > 0) {
      all_insertions[[gene_symbol]] <- result
      cat(" OK (", nrow(result), " insertions)\n")
    } else {
      failed_genes <- c(failed_genes, gene_symbol)
      cat(" FAILED\n")
    }
  }

  # Combine all results into one comprehensive table
  cat("\n=== Building Comprehensive Table ===\n")
  master_table <- do.call(rbind, all_insertions)

  # Reorder columns for clarity (now includes NEPC and RNA-seq info)
  master_table <- master_table[, c(
    "gene_name", "entrez_id", "gene_chromosome", "gene_strand", "gene_start", "gene_end",
    "cis_start", "cis_end",
    "sample", "is_nepc", "has_rnaseq", "rnaseq_id",
    "insertion_chr", "insertion_start", "insertion_end", "insertion_midpoint",
    "read_count", "insertion_strand", "distance_to_gene", "is_distant",
    "is_exonic", "is_same_strand", "classification"
  )]

  # Ensure tables directory exists
  dir.create(tables.dir, showWarnings = FALSE, recursive = TRUE)

  # Count genes with/without lollipop plots
  genes_with_nearby <- unique(master_table$gene_name[!master_table$is_distant])
  genes_all_distant <- setdiff(unique(master_table$gene_name), genes_with_nearby)

  cat("\n========================================\n")
  cat("ANALYSIS COMPLETE\n")
  cat("========================================\n")
  cat("Total insertions:", nrow(master_table), "\n")
  cat("Total genes:", length(unique(master_table$gene_name)), "\n")
  cat("  - Genes with nearby insertions (have plots):", length(genes_with_nearby), "\n")
  cat("  - Genes with only distant insertions (no plots):", length(genes_all_distant), "\n")
  cat("Failed genes:", length(failed_genes), "\n")
  cat("Skipped (no Entrez):", length(skipped_genes), "\n")
  if(generate_plots) cat("\nPlots saved to:", plots.dir, "\n")

  # Print classification summary
  cat("\n--- Classification Summary ---\n")
  print(table(master_table$classification))

  # Print NEPC/RNA-seq summary
  cat("\n--- Sample Phenotype Summary ---\n")
  cat("Insertions from NEPC samples:", sum(master_table$is_nepc, na.rm = TRUE), "\n")
  cat("Insertions from non-NEPC samples:", sum(!master_table$is_nepc, na.rm = TRUE), "\n")
  cat("Insertions with matched RNA-seq:", sum(master_table$has_rnaseq, na.rm = TRUE), "\n")

  return(list(
    master_table = master_table,
    failed = failed_genes,
    skipped = skipped_genes,
    genes_with_plots = genes_with_nearby,
    genes_no_plots = genes_all_distant
  ))
}

#' Process a single gene - extract insertions and classify them
#' Uses CIS region from TAPDANCE to find all insertions, then classifies based on
#' proximity to gene body. Insertions >2kb from gene body are classified as "distant".
#' @param gene_symbol Gene symbol
#' @param gene_entrez Entrez gene ID
#' @param distance_threshold Max distance from gene body for exonic/intronic classification (default 2000bp)
#' @param generate_plot Whether to generate PDF plot (only for genes with non-distant insertions)
#' @return Data frame with insertion details
process_gene <- function(gene_symbol, gene_entrez, distance_threshold = 2000,
                         generate_plot = TRUE) {

  # Get exons for this gene
  gene_exons <- all.exonic_parts[all.gene_id %in% gene_entrez]
  if(length(gene_exons) == 0) return(NULL)

  # Get strand and chromosome
  strand <- names(sort(table(strand(gene_exons)), decreasing=TRUE)[1])
  chromosome <- names(sort(table(seqnames(gene_exons)), decreasing=TRUE)[1])

  # Calculate gene boundaries
  gene_start_coord <- min(start(gene_exons))
  gene_end_coord <- max(end(gene_exons))

  # Get CIS data for this gene (includes CIS region from TAPDANCE)
  CIS <- CIS_.001 %>% filter(gene_name == gene_symbol)
  if(nrow(CIS) == 0) return(NULL)

  CIS_samples <- trimws(unlist(strsplit(CIS$library_name, "::")))
  CIS_samples <- CIS_samples[CIS_samples != ""]

  # Use CIS region from TAPDANCE to find ALL insertions for this gene
  cis_chr <- CIS$CIS_chr[1]
  cis_start <- CIS$CIS_start[1]
  cis_end <- CIS$CIS_end[1]

  # Find all insertions in CIS region
  CIS_subset <- simplified_bed %>%
    filter(V1 == cis_chr) %>%
    filter(V3 > cis_start) %>%
    filter(V2 < cis_end) %>%
    mutate(midpoint = (V2 + V3) / 2) %>%
    filter(V4 %in% CIS_samples)

  if(nrow(CIS_subset) == 0) return(NULL)

  # Calculate distance to gene body for each insertion
  dist_to_gene <- pmin(
    abs(CIS_subset$midpoint - gene_start_coord),
    abs(CIS_subset$midpoint - gene_end_coord)
  )
  # If insertion is within gene body, distance is 0
  within_gene <- CIS_subset$midpoint >= gene_start_coord & CIS_subset$midpoint <= gene_end_coord
  dist_to_gene[within_gene] <- 0

  # Determine if insertion is "distant" (>threshold from gene body)
  is_distant <- dist_to_gene > distance_threshold

  # Create insertions GRanges for exonic classification
  insertions.gr <- GRanges(
    seqnames = rep(chromosome, nrow(CIS_subset)),
    ranges = IRanges(CIS_subset$midpoint, width = 1)
  )

  # Classify: exonic vs intronic (only meaningful for non-distant insertions)
  is_exonic <- countOverlaps(insertions.gr, gene_exons) > 0

  # Classify: same strand vs opposite strand
  insertion_strand <- CIS_subset$V6
  if(strand == "+") {
    is_same_strand <- insertion_strand == "+"
  } else {
    is_same_strand <- insertion_strand == "-"
  }

  # Build classification - distant insertions get "distant" classification
  classification <- case_when(
    is_distant ~ "distant",
    is_exonic & is_same_strand ~ "exonic_same",
    is_exonic & !is_same_strand ~ "exonic_opposite",
    !is_exonic & is_same_strand ~ "intronic_same",
    !is_exonic & !is_same_strand ~ "intronic_opposite"
  )

  # Look up sample metadata for NEPC status and RNA-seq availability
  sample_ids <- CIS_subset$V4
  sample_info <- sample_metadata[match(sample_ids, sample_metadata$sample_id), ]

  # Handle missing samples (use defaults: NA)
  is_nepc <- sample_info$is_nepc
  has_rnaseq <- sample_info$has_rnaseq
  rnaseq_id <- sample_info$rnaseq_id

  # Build result data frame
  result_df <- data.frame(
    gene_name = gene_symbol,
    entrez_id = gene_entrez,
    gene_chromosome = chromosome,
    gene_strand = strand,
    gene_start = gene_start_coord,
    gene_end = gene_end_coord,
    cis_start = cis_start,
    cis_end = cis_end,
    sample = CIS_subset$V4,
    insertion_chr = CIS_subset$V1,
    insertion_start = CIS_subset$V2,
    insertion_end = CIS_subset$V3,
    insertion_midpoint = CIS_subset$midpoint,
    read_count = CIS_subset$V5,
    insertion_strand = CIS_subset$V6,
    distance_to_gene = dist_to_gene,
    is_distant = is_distant,
    is_exonic = is_exonic,
    is_same_strand = is_same_strand,
    classification = classification,
    is_nepc = is_nepc,
    has_rnaseq = has_rnaseq,
    rnaseq_id = rnaseq_id,
    stringsAsFactors = FALSE
  )

  # Generate plot only if there are non-distant insertions
  has_nearby_insertions <- any(!is_distant)
  if(generate_plot && has_nearby_insertions) {
    # Define plot region around gene body
    plot_start <- gene_start_coord - distance_threshold
    plot_end <- gene_end_coord + distance_threshold

    # Filter to non-distant insertions for plotting
    nearby_idx <- !is_distant
    CIS_nearby <- CIS_subset[nearby_idx, ]

    tryCatch({
      generate_lollipop_plot(gene_symbol, gene_entrez, gene_exons, strand,
                             chromosome, plot_start, plot_end,
                             CIS_nearby, is_exonic[nearby_idx],
                             is_same_strand[nearby_idx], classification[nearby_idx])
    }, error = function(e) {
      # Silently skip plot errors
    })
  }

  return(result_df)
}

#' Generate lollipop plot for a gene (separated from data processing)
#' Now includes NEPC status (shape) and RNA-seq availability (border width)
#' @param gene_symbol Gene symbol
#' @param gene_entrez Entrez gene ID
#' @param gene_exons GRanges of exons
#' @param strand Gene strand
#' @param chromosome Chromosome
#' @param start_bp Plot start position
#' @param end_bp Plot end position
#' @param CIS_subset Data frame with insertion data (must include V4=sample, V5=score, midpoint)
#' @param is_exonic Logical vector - is insertion exonic
#' @param is_same_strand Logical vector - is insertion same strand as gene
#' @param classification Character vector of classifications
generate_lollipop_plot <- function(gene_symbol, gene_entrez, gene_exons, strand,
                                    chromosome, start_bp, end_bp,
                                    CIS_subset, is_exonic, is_same_strand, classification) {

  # Style exons
  gene_exons$fill <- COLORS$exon_fill
  gene_exons$color <- COLORS$exon_border
  gene_exons$height <- 0.08
  names(gene_exons) <- NULL

  # TSS marker
  if(strand == "+") {
    tss_pos <- min(start(gene_exons))
  } else {
    tss_pos <- max(end(gene_exons))
  }
  tss_marker <- GRanges(seqnames = chromosome,
                        ranges = IRanges(start = tss_pos, width = 50),
                        strand = strand)
  tss_marker$fill <- COLORS$promoter_fill
  tss_marker$color <- COLORS$promoter_border
  tss_marker$height <- 0.12
  names(tss_marker) <- NULL

  # Range for plotting
  range <- GRanges(seqnames = chromosome,
                   ranges = IRanges(start = start_bp, end = end_bp),
                   strand = strand)

  # Create insertions GRanges with styling
  insertions.gr <- GRanges(
    seqnames = rep(chromosome, nrow(CIS_subset)),
    ranges = IRanges(CIS_subset$midpoint, width = 1)
  )

  # Assign colors based on exonic/intronic and strand

insertions.gr$color <- case_when(
    is_exonic & is_same_strand ~ COLORS$exonic_same,
    is_exonic & !is_same_strand ~ COLORS$exonic_opposite,
    !is_exonic & is_same_strand ~ COLORS$intronic_same,
    !is_exonic & !is_same_strand ~ COLORS$intronic_opposite
  )

  # Look up sample metadata for NEPC status and RNA-seq availability
  sample_ids <- CIS_subset$V4
  sample_info <- sample_metadata[match(sample_ids, sample_metadata$sample_id), ]

  # Handle missing samples (use defaults)
  is_nepc <- ifelse(is.na(sample_info$is_nepc), FALSE, sample_info$is_nepc)
  has_rnaseq <- ifelse(is.na(sample_info$has_rnaseq), FALSE, sample_info$has_rnaseq)

  # Shape: circle (pch 21) = NEPC, square (pch 22) = non-NEPC
  # trackViewer uses shape names: "circle", "square", "diamond", "triangle_point_up", etc.
  insertions.gr$shape <- ifelse(is_nepc, "circle", "square")

  # Border width: thick (lwd 2.5) = has RNA-seq, thin (lwd 0.8) = no RNA-seq
  insertions.gr$lwd <- ifelse(has_rnaseq, 2.5, 0.8)

  insertions.gr$score <- CIS_subset$V5
  insertions.gr$height <- CIS_subset$V5
  insertions.gr$border <- "#2c3e50"
  insertions.gr$alpha <- 0.85
  insertions.gr$cex <- 0.9  # Slightly larger for visibility
  insertions.gr$SNPsideID <- "top"

  # X-axis labels
  xaxis <- floor(seq(start_bp, end_bp, length.out = 6))
  names(xaxis) <- format(xaxis, big.mark = ",", scientific = FALSE)

  # Ensure output directory exists
  dir.create(plots.dir, showWarnings = FALSE, recursive = TRUE)

  # Generate plot
  output_file <- file.path(plots.dir, paste0(gene_symbol, ".pdf"))
  pdf(output_file, width = 10, height = 5)  # Slightly taller for two-row legend

  features <- c(gene_exons, tss_marker)

  lolliplot(insertions.gr, features, range,
            ylab = "Read Count",
            ylab.gp = gpar(fontsize = 11, fontface = "bold"),
            yaxis.gp = gpar(fontsize = 9),
            xaxis.gp = gpar(fontsize = 8),
            xaxis = xaxis,
            legend = FALSE,
            jitter = "node")

  # Custom legend - Row 1: Colors (exonic/intronic status)
  pushViewport(viewport(x = 0.5, y = 0.88, width = 0.9, height = 0.04))
  grid.draw(legendGrob(
    labels = c("Exonic (sense)", "Exonic (antisense)",
               "Intronic (sense)", "Intronic (antisense)"),
    pch = 21,
    gp = gpar(col = c(COLORS$exonic_same, COLORS$exonic_opposite,
                      COLORS$intronic_same, COLORS$intronic_opposite),
              fill = c(COLORS$exonic_same, COLORS$exonic_opposite,
                       COLORS$intronic_same, COLORS$intronic_opposite)),
    ncol = 4, byrow = TRUE
  ))
  popViewport()

  # Custom legend - Row 2: Shapes (NEPC status) and borders (RNA-seq)
  pushViewport(viewport(x = 0.5, y = 0.82, width = 0.9, height = 0.04))
  grid.draw(legendGrob(
    labels = c("NEPC", "non-NEPC", "Has RNA-seq", "No RNA-seq"),
    pch = c(21, 22, 21, 21),
    gp = gpar(col = c("gray40", "gray40", "black", "gray60"),
              fill = c("gray70", "gray70", "gray70", "gray70"),
              lwd = c(1.5, 1.5, 2.5, 0.8)),
    ncol = 4, byrow = TRUE
  ))
  popViewport()

  # Title
  grid.text(paste0(gene_symbol, " (", strand, " strand)"),
            x = 0.5, y = 0.94,
            gp = gpar(cex = 1.4, fontface = "bold.italic"))

  dev.off()
}

#' Compare our classification with Sai's analysis
#' @param master_table Output from run_full_analysis()$master_table
#' @param sai_file Path to Sai's insertions_with_type.csv
compare_with_sai <- function(master_table,
                              sai_file = "/Users/matteodibernardo/Library/CloudStorage/Dropbox/Matteo_revisions/Sai_Analysis_from_email/insertions_with_type.csv") {

  cat("\n=== Comparing with Sai's Classification ===\n")

  # Load Sai's data
  sai_data <- read.csv(sai_file)
  cat("Loaded", nrow(sai_data), "insertions from Sai's file\n")
  cat("Our data has", nrow(master_table), "insertions\n")

  # Create matching key (chr:start:end:sample)
  master_table$match_key <- paste(
    master_table$insertion_chr,
    master_table$insertion_start,
    master_table$insertion_end,
    master_table$sample,
    sep = ":"
  )

  sai_data$match_key <- paste(
    sai_data$Chromosome,
    sai_data$Start,
    sai_data$End,
    sai_data$sample,
    sep = ":"
  )

  # Merge
  merged <- merge(
    master_table[, c("match_key", "gene_name", "is_exonic", "classification")],
    sai_data[, c("match_key", "insertion_type")],
    by = "match_key",
    all = FALSE
  )

  cat("Matched", nrow(merged), "insertions\n")

  # Compare: our exonic should match Sai's "exonic"
  merged$our_exonic <- merged$is_exonic
  merged$sai_exonic <- merged$insertion_type == "exonic"
  merged$match <- merged$our_exonic == merged$sai_exonic

  n_match <- sum(merged$match)
  n_mismatch <- sum(!merged$match)
  agreement_rate <- round(n_match / nrow(merged) * 100, 2)

  cat("\n--- Classification Agreement ---\n")
  cat("Matching:", n_match, "\n")
  cat("Mismatching:", n_mismatch, "\n")
  cat("Agreement rate:", agreement_rate, "%\n")

  # Save comparison
  dir.create(tables.dir, showWarnings = FALSE, recursive = TRUE)
  output_file <- file.path(tables.dir, "classification_comparison.csv")
  write.csv(merged, output_file, row.names = FALSE)
  cat("\nComparison saved to:", output_file, "\n")

  if(n_mismatch > 0 && n_mismatch <= 20) {
    cat("\n--- Mismatches ---\n")
    print(merged[!merged$match, c("gene_name", "match_key", "our_exonic", "sai_exonic")])
  }

  return(merged)
}

# ============================================================================
# SLEEPING BEAUTY MECHANISM CLASSIFICATION
# ============================================================================
# Based on SB transposon biology:
# - Promoter-proximal + sense strand = ACTIVATION (MSCV promoter drives transcription)
# - Early intronic + sense strand = ACTIVATION (splice into downstream exons)
# - Everything else = DISRUPTION (truncation, antisense, etc.)

#' Classify insertions based on Sleeping Beauty mechanism
#'
#' Rules for ACTIVATION-capable insertions:
#' 1. Promoter-proximal: -2000bp to +500bp from TSS AND sense strand
#' 2. Early intronic: In intron 1-3 OR first 20% of transcript AND sense strand
#' All other insertions are classified as DISRUPTION.
#'
#' @param master_table Output from run_full_analysis()$master_table
#' @return Updated master_table with sb_mechanism column and summary
classify_sb_mechanism <- function(master_table = NULL) {

  cat("\n========================================\n")
  cat("SLEEPING BEAUTY MECHANISM CLASSIFICATION\n")
  cat("========================================\n")

  # Load master table if not provided
  if(is.null(master_table)) {
    stop("Master table must be provided. Run run_full_analysis() first.")
  }

  # Ensure data is loaded for exon queries
  if(!exists("all.exonic_parts")) {
    load_lollipop_data()
  }

  cat("\nClassification rules:\n")
  cat("  ACTIVATION: Promoter-proximal (-2kb to +500bp TSS, sense) OR early-intronic (intron 1-3 or first 20%, sense)\n")
  cat("  DISRUPTION: Everything else\n\n")

  # Initialize columns
  master_table$tss_position <- NA
  master_table$distance_from_tss <- NA
  master_table$is_promoter_proximal <- FALSE
  master_table$intron_number <- NA
  master_table$pct_into_transcript <- NA
  master_table$is_early_intronic <- FALSE
  master_table$sb_mechanism <- "disruption"

  # Process each gene
  genes <- unique(master_table$gene_name)
  cat("Processing", length(genes), "genes...\n")

  for(gene in genes) {
    gene_idx <- master_table$gene_name == gene
    gene_data <- master_table[gene_idx, ]
    gene_entrez <- gene_data$entrez_id[1]
    gene_strand <- gene_data$gene_strand[1]
    gene_start <- gene_data$gene_start[1]
    gene_end <- gene_data$gene_end[1]

    # Get exons for this gene
    gene_exons <- all.exonic_parts[all.gene_id %in% gene_entrez]
    if(length(gene_exons) == 0) next

    # Calculate TSS position
    if(gene_strand == "+") {
      tss_pos <- min(start(gene_exons))
      transcript_end <- max(end(gene_exons))
    } else {
      tss_pos <- max(end(gene_exons))
      transcript_end <- min(start(gene_exons))
    }

    transcript_length <- abs(transcript_end - tss_pos)

    # Store TSS position
    master_table$tss_position[gene_idx] <- tss_pos

    # Get exon coordinates sorted by position in transcript
    exon_df <- as.data.frame(gene_exons)
    if(gene_strand == "+") {
      exon_df <- exon_df[order(exon_df$start), ]
    } else {
      exon_df <- exon_df[order(exon_df$end, decreasing = TRUE), ]
    }

    # Calculate intron boundaries
    # Intron N is between exon N and exon N+1
    n_exons <- nrow(exon_df)
    intron_ranges <- list()
    if(n_exons > 1) {
      for(i in 1:(n_exons - 1)) {
        if(gene_strand == "+") {
          intron_start <- exon_df$end[i] + 1
          intron_end <- exon_df$start[i + 1] - 1
        } else {
          intron_start <- exon_df$end[i + 1] + 1
          intron_end <- exon_df$start[i] - 1
        }
        if(intron_end >= intron_start) {
          intron_ranges[[i]] <- c(start = intron_start, end = intron_end)
        }
      }
    }

    # Classify each insertion for this gene
    for(i in which(gene_idx)) {
      ins_mid <- master_table$insertion_midpoint[i]
      ins_strand <- master_table$insertion_strand[i]
      is_same_strand <- master_table$is_same_strand[i]
      is_exonic <- master_table$is_exonic[i]
      is_distant <- master_table$is_distant[i]

      # Distance from TSS (signed: negative = upstream, positive = downstream)
      if(gene_strand == "+") {
        dist_from_tss <- ins_mid - tss_pos
      } else {
        dist_from_tss <- tss_pos - ins_mid  # Reverse for minus strand
      }
      master_table$distance_from_tss[i] <- dist_from_tss

      # Percent into transcript (0% = TSS, 100% = end)
      # Keep raw value (can be negative for upstream, >100 for downstream)
      if(gene_strand == "+") {
        pct_into <- (ins_mid - tss_pos) / transcript_length * 100
      } else {
        pct_into <- (tss_pos - ins_mid) / transcript_length * 100
      }
      master_table$pct_into_transcript[i] <- round(pct_into, 2)

      # DISTANT insertions are automatically DISRUPTION - skip further classification
      if(is_distant) {
        master_table$is_promoter_proximal[i] <- FALSE
        master_table$is_early_intronic[i] <- FALSE
        master_table$sb_mechanism[i] <- "disruption"
        next
      }

      # Check promoter-proximal: -2000 to +500 bp from TSS AND sense
      is_promoter_proximal <- (dist_from_tss >= -2000 & dist_from_tss <= 500) & is_same_strand
      master_table$is_promoter_proximal[i] <- is_promoter_proximal

      # Determine intron number (if intronic)
      intron_num <- NA
      if(!is_exonic && length(intron_ranges) > 0) {
        for(j in seq_along(intron_ranges)) {
          if(!is.null(intron_ranges[[j]])) {
            if(ins_mid >= intron_ranges[[j]]["start"] && ins_mid <= intron_ranges[[j]]["end"]) {
              intron_num <- j
              break
            }
          }
        }
      }
      master_table$intron_number[i] <- intron_num

      # Check early-intronic: intron 1-3 OR first 20% of transcript AND sense
      # Must be WITHIN the gene body (pct_into between 0-100)
      is_early_intronic <- FALSE
      if(!is_exonic && is_same_strand && pct_into >= 0 && pct_into <= 100) {
        in_early_intron <- !is.na(intron_num) && intron_num <= 3
        in_first_20pct <- pct_into <= 20
        is_early_intronic <- in_early_intron | in_first_20pct
      }
      master_table$is_early_intronic[i] <- is_early_intronic

      # Final classification
      if(is_promoter_proximal | is_early_intronic) {
        master_table$sb_mechanism[i] <- "activation"
      } else {
        master_table$sb_mechanism[i] <- "disruption"
      }
    }
  }

  # Generate summary
  cat("\n=== SB MECHANISM CLASSIFICATION RESULTS ===\n\n")

  # Overall counts
  cat("Overall Insertion Classification:\n")
  print(table(master_table$sb_mechanism))
  cat("\n")

  # By original classification
  cat("SB Mechanism by Original Classification:\n")
  print(table(master_table$classification, master_table$sb_mechanism))
  cat("\n")

  # Breakdown of activation insertions
  activation_ins <- master_table[master_table$sb_mechanism == "activation", ]
  cat("Activation-capable insertions breakdown:\n")
  cat("  Promoter-proximal:", sum(activation_ins$is_promoter_proximal), "\n")
  cat("  Early-intronic:", sum(activation_ins$is_early_intronic & !activation_ins$is_promoter_proximal), "\n")
  cat("  Both:", sum(activation_ins$is_promoter_proximal & activation_ins$is_early_intronic), "\n")
  cat("\n")

  # Gene-level summary
  gene_summary <- master_table %>%
    group_by(gene_name) %>%
    summarize(
      n_insertions = n(),
      n_activation = sum(sb_mechanism == "activation"),
      n_disruption = sum(sb_mechanism == "disruption"),
      n_distant = sum(classification == "distant"),
      n_exonic = sum(is_exonic & classification != "distant"),
      has_activation = any(sb_mechanism == "activation"),
      has_exonic = any(is_exonic & classification != "distant"),
      .groups = "drop"
    )

  cat("Gene-level Summary:\n")
  cat("  Total genes:", nrow(gene_summary), "\n")
  cat("  Genes with activation-capable insertions:", sum(gene_summary$has_activation), "\n")
  cat("  Genes with exonic insertions:", sum(gene_summary$has_exonic), "\n")
  cat("  Genes with only distant insertions:", sum(gene_summary$n_insertions == gene_summary$n_distant), "\n")
  cat("\n")

  # Ensure tables directory exists
  dir.create(tables.dir, showWarnings = FALSE, recursive = TRUE)

  # Save master table (all insertions)
  output_file <- file.path(tables.dir, "all_insertions_master.csv")
  write.csv(master_table, output_file, row.names = FALSE)
  cat("Master table saved to:", output_file, "\n")

  # Save gene summary
  summary_file <- file.path(tables.dir, "gene_summary.csv")
  write.csv(gene_summary, summary_file, row.names = FALSE)
  cat("Gene summary saved to:", summary_file, "\n")

  # Save key genes subset (reviewer-requested)
  key_genes_table <- master_table[master_table$gene_name %in% KEY_GENES, ]
  key_genes_file <- file.path(tables.dir, "key_genes_insertions.csv")
  write.csv(key_genes_table, key_genes_file, row.names = FALSE)
  cat("Key genes table saved to:", key_genes_file, "(", nrow(key_genes_table), "insertions across",
      length(unique(key_genes_table$gene_name)), "genes)\n")

  return(list(
    master_table = master_table,
    gene_summary = gene_summary,
    key_genes_table = key_genes_table
  ))
}

#' Generate overall summary statistics table for the paper
#' @param master_table Output from classify_sb_mechanism()$master_table
#' @return Summary table suitable for paper
generate_paper_summary <- function(master_table = NULL) {

  cat("\n========================================\n")
  cat("GENERATING PAPER SUMMARY TABLE\n")
  cat("========================================\n")

  # Load master table if not provided
  if(is.null(master_table)) {
    master_file <- file.path(tables.dir, "all_insertions_master.csv")
    if(!file.exists(master_file)) {
      stop("Master table not found. Run run_full_analysis() and classify_sb_mechanism() first.")
    }
    master_table <- read.csv(master_file, stringsAsFactors = FALSE)
  }

  # Check if sb_mechanism exists
  has_sb <- "sb_mechanism" %in% colnames(master_table)

  # Calculate summary statistics
  total_cis_genes <- length(unique(master_table$gene_name))
  total_insertions <- nrow(master_table)

  # Genes with nearby vs distant only
  genes_nearby <- unique(master_table$gene_name[master_table$classification != "distant"])
  genes_distant_only <- setdiff(unique(master_table$gene_name), genes_nearby)

  # Insertion location breakdown (excluding distant)
  nearby_ins <- master_table[master_table$classification != "distant", ]
  n_exonic <- sum(nearby_ins$is_exonic)
  n_intronic <- sum(!nearby_ins$is_exonic)

  # Strand orientation
  n_same_strand <- sum(nearby_ins$is_same_strand)
  n_opposite_strand <- sum(!nearby_ins$is_same_strand)

  # Build summary table
  summary_df <- data.frame(
    Category = c(
      "Total CIS genes",
      "Total insertions",
      "Genes with gene-body insertions (≤2kb)",
      "Genes with only distant insertions (>2kb)",
      "",
      "Gene-body insertions:",
      "  Exonic insertions",
      "  Intronic insertions",
      "  Same-strand insertions (sense)",
      "  Opposite-strand insertions (antisense)"
    ),
    Count = c(
      total_cis_genes,
      total_insertions,
      length(genes_nearby),
      length(genes_distant_only),
      NA,
      nrow(nearby_ins),
      n_exonic,
      n_intronic,
      n_same_strand,
      n_opposite_strand
    ),
    stringsAsFactors = FALSE
  )

  # Add SB mechanism if available
  if(has_sb) {
    sb_rows <- data.frame(
      Category = c(
        "",
        "SB Mechanism Classification:",
        "  Activation-capable insertions",
        "  Disruption insertions",
        "",
        "Genes with activation-capable insertions",
        "Genes with exonic insertions"
      ),
      Count = c(
        NA,
        NA,
        sum(master_table$sb_mechanism == "activation"),
        sum(master_table$sb_mechanism == "disruption"),
        NA,
        length(unique(master_table$gene_name[master_table$sb_mechanism == "activation"])),
        length(unique(nearby_ins$gene_name[nearby_ins$is_exonic]))
      ),
      stringsAsFactors = FALSE
    )
    summary_df <- rbind(summary_df, sb_rows)
  }

  # Print to console
  cat("\n")
  print(summary_df, row.names = FALSE)

  # Save
  dir.create(tables.dir, showWarnings = FALSE, recursive = TRUE)
  output_file <- file.path(tables.dir, "paper_summary.csv")
  write.csv(summary_df, output_file, row.names = FALSE)
  cat("\nSummary saved to:", output_file, "\n")

  return(summary_df)
}

# ============================================================================
# USAGE
# ============================================================================
# Full workflow:
#   source("libs/lollipop.R")
#   results <- run_full_analysis()                    # Generate plots and raw table
#   sb_results <- classify_sb_mechanism(results$master_table)  # Add SB mechanism, save tables
#   summary <- generate_paper_summary(sb_results$master_table) # Paper summary stats
#
# Optional: Validate against Sai's classification (run once, then remove call)
#   comparison <- compare_with_sai(results$master_table)
#
# Output structure:
#   experiments/lollipop/
#   ├── plots/                     # PDF lollipop plots (172 files)
#   │   ├── Sirt1.pdf
#   │   └── ...
#   └── tables/                    # CSV data tables
#       ├── all_insertions_master.csv    # All 3079 insertions with SB mechanism
#       ├── gene_summary.csv             # Per-gene summary (282 genes)
#       ├── key_genes_insertions.csv     # 15 reviewer-requested genes only
#       └── paper_summary.csv            # Overall statistics for paper
#
# To run for a single gene (with plot):
#   load_lollipop_data()
#   result <- lolliplot_updated("Sirt1", "93759", promoter_bp = 2000)
