reports.dir <- "experiments/lollipop"

library(trackViewer)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(dplyr)
library(biomaRt)     

# Read in the data
CIS_.001 <- read.csv("experiments/table/TAPDANCE_001.csv")
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
all.exonic_parts <- exonicParts(txdb, linked.to.single.gene.only=TRUE)
all.gene_id <- mcols(all.exonic_parts)$gene_id
all.genes <- genes(txdb)
all.promoters <- promoters(all.genes,upstream=100,downstream=50)
simplified_bed <- read.table("TAPDANCE_NEPC/results/all/all-nr-prostate-0.001.BED", sep = "\t",header = FALSE, skip = 1)

# Function to create a lollipop plot
lolliplot_custom <- function(gene_symbol, gene_entrez, promoter_bp = 5000, overhang_bp = 0, names = TRUE){
  # Print gene information
  cat("Generating lollipop for gene symbol:", gene_symbol, "entrez ID:", gene_entrez, "\n")
  
  # Get the exons
  gene_exons <- all.exonic_parts[all.gene_id %in% gene_entrez]
  cat("Number of exons found:", length(gene_exons), "\n")
  gene_exons$fill <- "grey"
  gene_exons$height <- .02
  gene_exons$color<- "grey"
  names(gene_exons) <- rep("Exons", length(gene_exons))  # Add this line
  
  # Get the promoter
  promoter <- all.promoters[all.promoters$gene_id %in% gene_entrez]
  cat("Number of promoters found:", length(promoter), "\n")
  promoter$fill <- "red"
  promoter$height <- .05
  promoter$color<- "red"
  names(promoter) <- "Promoter"  # Add this line
  
  # Get the strand and chromosome number
  strand <- names(sort(table(strand(gene_exons)),decreasing=TRUE)[1])
  chromosome_number <- names(sort(table(seqnames(gene_exons)),decreasing=TRUE)[1])
  
  # Set the gene start and end based on strand
  if(strand == "+"){
    gene_start <- min(start(gene_exons)) - promoter_bp
    gene_end <- max(end(gene_exons)) + overhang_bp
  }
  if(strand == "-"){
    gene_start <- max(end(gene_exons)) + promoter_bp
    gene_end <- min(start(gene_exons)) - overhang_bp
  }
  start_bp <- min(gene_start,gene_end)
  end_bp <- max(gene_start,gene_end)

  # Set the range based on the gene start and end
  range <-GRanges(seqnames=chromosome_number,
                  ranges=IRanges(start = start_bp-promoter_bp, end = end_bp+promoter_bp),
                  strand=strand)

  # Extract the relevant insertions
  CIS <- CIS_.001 %>% dplyr::filter(gene_name == gene_symbol)
  cat("CIS entries found:", nrow(CIS), "\n")
  CIS_samples <- trimws(unlist(strsplit(CIS$library_name, "::")), which = "both")
  CIS_samples <- CIS_samples[CIS_samples != ""]

  # Subset the insertions
  CIS_subset <- simplified_bed %>% dplyr::filter(V1 == chromosome_number) %>%
    dplyr::filter(V3 > start_bp) %>%
    dplyr::filter(V2 < end_bp) %>%
    dplyr::mutate(midpoint = (V2+V3)/2) %>%
    dplyr::filter(V4 %in% CIS_samples)

  if(nrow(CIS_subset) == 0) {
    cat("\nWARNING: No insertions found! Debugging info:")
    cat("\n - Looking for samples:", paste(CIS_samples, collapse=", "))
    cat("\n - In chromosome:", chromosome_number)
    cat("\n - Between positions:", start_bp, "and", end_bp, "\n")
    return(NULL)
  }
  
  # Construct the insertions GRanges object
  insertions <- CIS_subset$midpoint
  if(names){
    insertions.gr <- GRanges(seqnames = rep(chromosome_number, length(insertions)), 
                             ranges = IRanges(insertions, width=1), 
                             height = CIS_subset$V5,
                             names = CIS_subset$V4)
  } else{
    insertions.gr <- GRanges(seqnames = rep(chromosome_number, length(insertions)), 
                             ranges = IRanges(insertions, width=1),
                             height = CIS_subset$V5)
  }
  insertions.gr$score <- insertions.gr$height
  insertions.gr$SNPsideID <- ifelse(CIS_subset$V6 == "+","top","top")
  insertions.gr$cex <- .3
  insertions.gr$alpha <- .9
  insertions.gr$label.parameter.gp <- gpar(col = "grey")
  insertions.gr$label.parameter.gp <- gpar(cex = .5)
  if(strand == "+"){
    insertions.gr$color <- ifelse(CIS_subset$V6 == "+","green","red")
    insertions.gr$shape <- ifelse(CIS_subset$V6 == "+","triangle_point_down","triangle_point_up")
  }
  if(strand == "-"){
    insertions.gr$color <- ifelse(CIS_subset$V6 == "-","green","red")
    insertions.gr$shape <- ifelse(CIS_subset$V6 == "-","triangle_point_down","triangle_point_up")
  }

  # Lolliplot labels
  xaxis <- floor(seq(start_bp-100000/30, end_bp+100000/30, length.out = 10)[2:9]) ## define the position
  names(xaxis) <- xaxis 
  
  # Create the lollipop plot
  pdf(file.path(reports.dir, paste(gene_symbol,".pdf",sep = "")), height=5)
  lolliplot(insertions.gr, c(gene_exons,promoter), range,
            ylab = "transposon read count",
            yaxis.gp = gpar(fontsize = 6),
            xaxis.gp = gpar(fontsize = 4.5),
            ylab.gp = gpar(fontsize = 10),
            xaxis = xaxis,
  )
  # Adjust title position
  grid.text(gene_symbol, just="top", y = .85,  # moved down from .65
            gp=gpar(cex=1.5, fontface="italic"))
  
  # Create viewport for legend
  pushViewport(viewport(x = 0.2, y = 0.2, just = "left"))  # moved to left side to avoid overlapping with tall peaks
  
  grid.polygon(x = c(0.046, 0.05, 0.054),  
               y = c(0.95, 0.94, 0.95), 
               gp = gpar(fill = "green", col = "black"))  # added black border
  
  # Red upward triangle with black border           
  grid.polygon(x = c(0.046, 0.05, 0.054),  
               y = c(0.90, 0.91, 0.90), 
               gp = gpar(fill = "red", col = "black"))  # added black border
  
  # Adjust text position to align with triangles
  grid.text("Same strand", x = 0.06, y = 0.945,    # aligned with green triangle
            just = "left", gp = gpar(fontsize = 6))
  grid.text("Opposite strand", x = 0.06, y = 0.905, # aligned with red triangle
            just = "left", gp = gpar(fontsize = 6))
  
  popViewport()
  
  dev.off()

  #write a csv for the genes relevant insertions
  CIS_subset %>% write.csv(file.path(reports.dir,paste(gene_symbol,".csv",sep = "")))
}

# Get the gene IDs as Entrez
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host="useast.ensembl.org")
mouse_gene_ids  <- CIS_.001$gene_name
conversion <- getBM(attributes=c('mgi_symbol',
                                 'entrezgene_id'),
                    filters = 'mgi_symbol',
                    values = mouse_gene_ids,
                    mart = ensembl)
conversion <- conversion[!is.na(conversion$entrezgene_id),]

# Create the lollipop plots (impossible to eliminate the block legend, better to do this via photoshop)
for(i in 1:nrow(conversion)){
  tryCatch({
    lolliplot_custom(conversion$mgi_symbol[i], conversion$entrezgene_id[i], promoter_bp = 10000)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

# Some genes of interest may need increased promoter regions or an extra base pairs at the end to catch the cis
lolliplot_custom("Akap8","56399",promoter_bp = 10000,overhang_bp = 150000)
lolliplot_custom("Cxxc5","67393",promoter_bp = 50000)
lolliplot_custom("Evi2a","14017",promoter_bp = 50000)
