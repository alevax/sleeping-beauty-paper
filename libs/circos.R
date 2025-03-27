reports.dir <- "experiments/circos"

library(stringr)
library(biomaRt)
library(dplyr)
library(circlize)
library(ComplexHeatmap)
library(RColorBrewer)

# Read in the TAPDANCE output and get the gene annotation
TD <- read.csv("experiments/table/TAPDANCE_001.csv")

# Get the gene annotation
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl", host="https://useast.ensembl.org")

# Get the gene annotation
TD_all <- getBM(attributes = c("mgi_symbol","chromosome_name",'strand','transcript_start','transcript_end'), 
                mart = mouse,
                filters = "mgi_symbol",
                values = TD$gene_name,unique = TRUE)

# Get the unique gene annotation
TD_unique <- TD_all %>% 
  dplyr::mutate(length = transcript_end - transcript_start) %>% 
  dplyr::group_by(mgi_symbol,chromosome_name,strand) %>%
  dplyr::slice(which.max(length)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(chromosome_name %in% c(1:19,"X")) %>%
  dplyr::select(mgi_symbol,chromosome_name,transcript_start,transcript_end)

# Get the genes that are not in the unique gene annotation
left_out_TD <- TD %>% dplyr::filter(!(TD$gene_name %in% TD_unique$mgi_symbol))
for(i in 1:nrow(left_out_TD)){
  left_out_TD$chromosome_name[i] <- str_sub(unlist(strsplit(left_out_TD$pos[i], ":"))[1],4)
  left_out_TD$transcript_start[i] <- as.numeric(
    unlist(strsplit(unlist(strsplit(left_out_TD$pos[i], ":"))[2], "-"))[1])
  left_out_TD$transcript_end[i] <- as.numeric(
    unlist(strsplit(unlist(strsplit(left_out_TD$pos[i], ":"))[2], "-"))[2])
}

left_out_TD <- left_out_TD %>% dplyr::select(mgi_symbol = gene_name,chromosome_name,transcript_start,transcript_end)

# Combine the unique and left out gene annotation
TD_annotation_final <- rbind(left_out_TD,TD_unique)

# Read in the BED file
final_bed <- read.table("data/tapdance/raw_prostate.BED", sep = "\t",header = FALSE, skip = 1)
positive_bed <- final_bed %>% dplyr::filter(V6 == "+") %>% dplyr::select(V1,V2,V3,V5)
negative_bed <- final_bed %>% dplyr::filter(V6 == "-") %>% dplyr::select(V1,V2,V3,V5)

# Combine the BED files
bed_list <- list(positive_bed, negative_bed)  

# Get the unique gene annotation
TD_annotation_final <- TD_annotation_final %>% dplyr::mutate(chromosome_name = paste("chr",chromosome_name,sep="")) %>%
  dplyr::select(chromosome_name,transcript_start,transcript_end,mgi_symbol) 

for(i in 1:nrow(TD)){
  TD$chromosome_name[i] <- str_sub(unlist(strsplit(TD$pos[i], ":"))[1],4)
  TD$transcript_start[i] <- as.numeric(
    unlist(strsplit(unlist(strsplit(TD$pos[i], ":"))[2], "-"))[1])
  TD$transcript_end[i] <- as.numeric(
    unlist(strsplit(unlist(strsplit(TD$pos[i], ":"))[2], "-"))[2])
}

bed_TD <- TD %>% dplyr::mutate(chromosome_name = paste("chr",chromosome_name,sep="")) %>%
  dplyr::select(chr = chromosome_name,start = transcript_start,end = transcript_end) %>% unique()
bed_TD$value1 <- runif(n = nrow(bed_TD), min = -1, max = 1)
bed_TD$start <- bed_TD$start - 350000
bed_TD$end <- bed_TD$end + 350000
col_fun1 = colorRamp2(breaks = c(-1, 1), colors = c("blue", "red"))

# Read in the co-occurrence data
SB_co <- read.csv("experiments/table/SB_co-occurrence_CIS_genes.csv")[,2:4] %>% dplyr::filter(pValue < .01) %>% unique()

SB_co_1 <- SB_co %>% dplyr::select(gene_name = gene1)
SB_co_2 <- SB_co %>% dplyr::select(gene_name = gene2)

co_TD1 <- SB_co_1 %>% left_join(TD) %>% dplyr::mutate(chromosome_name = paste("chr",chromosome_name,sep="")) %>%
  dplyr::select(chr = chromosome_name,start = transcript_start,end = transcript_end) 
co_TD1$start <- co_TD1$start - 1000000
co_TD1$end <- co_TD1$end + 1000000

co_TD2 <- SB_co_2 %>% left_join(TD) %>% dplyr::mutate(chromosome_name = paste("chr",chromosome_name,sep="")) %>%
  dplyr::select(chr = chromosome_name,start = transcript_start,end = transcript_end) 
co_TD2$start <- co_TD2$start - 1000000
co_TD2$end <- co_TD2$end + 1000000

# Remove the co-occurrences that are not in the unique gene annotation
save <- vector()
for(i in 1:nrow(co_TD1)){
  save[i] <- ifelse(identical(co_TD1[i,], co_TD2[i,]),FALSE,TRUE)
}

co_TD1 <- co_TD1[which(save),]
co_TD2 <- co_TD2[which(save),]

pval.breaks <- seq(0, .01, length.out=100)
colrampcustom <- colorRampPalette(c("red", "white"))
pval.cols <- colrampcustom(100)
col_fun2 = colorRamp2(breaks = pval.breaks, colors = pval.cols)

circos.clear()

# Create the circos plot
pdf(file.path(reports.dir, "circos.pdf"), width = 47, height = 47)
# Set the circos parameters
circos.par(circle.margin= c(0.01,0.01,0.01,0.01), gap.after = c(rep(1, 20), 20),start.degree = 350)
# Initialize the circos plot
circos.initializeWithIdeogram(species = "mm10", 
                              plotType = NULL)
# Add the genomic labels
circos.genomicLabels(TD_annotation_final,
                     labels.column = 4, cex = 2.25,
                     side = "outside",
                     line_col = "grey",
                     connection_height = mm_h(100))
# Add the genomic ideogram
circos.genomicIdeogram(species = "mm10",track.height = mm_h(5))
# Add the genomic track
circos.track(track.index =3, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2], 
              gsub("chr", "", CELL_META$sector.index), adj = c(0.5, -0.8),cex = 4)
})
# Add the insertion positions
circos.genomicTrack(bed_TD, stack = TRUE, track.height = .1, bg.lwd = .5, panel.fun = function(region, value, ...) {
  circos.genomicRect(region, value, col = "black", border = NA, ...)
})

# Add the positive and negative strands
circos.genomicTrack(positive_bed,track.height = .1,bg.lwd = .5,panel.fun = function(region, value, ...) {
  circos.genomicPoints(region, value,pch = 20, cex = 1, col = "red",...)
})
circos.yaxis("left", at = 50000, labels.cex = 1.5, labels.niceFacing = F)

circos.genomicTrack(negative_bed,track.height = .1,bg.lwd = .5,panel.fun = function(region, value, ...) {
  circos.genomicPoints(region, value,pch = 20, cex = 1, col = "blue",...)
})
# Add the co-occurrences
circos.yaxis("left", at = 50000, labels.cex = 1.5, labels.niceFacing = F)
circos.genomicLink(co_TD1, co_TD2, col = col_fun2(SB_co$pValue[which(save)]), lwd = 5, 
                   border = NA)
# Add the legend
lgd_links = Legend(at = c(1e-25, 1e-3, .01), 
                   col_fun = col_fun2, 
                   title_position = "topcenter", 
                   labels_gp = gpar(fontsize = 15),
                   grid_height = unit(10, "mm"), grid_width = unit(10, "mm"),
                   title = "CIS \n co-occurrence \n p-value", 
                   title_gp = gpar(fontsize = 20, fontface = "bold.italic"), 
                   title_gap = unit(5.5, "mm"),
                   border = "black", background = "white")
lgd_points = Legend(at = c("negative", "positive"), 
                    type = "points", 
                    legend_gp = gpar(col = c("blue","red")), 
                    labels_gp = gpar(fontsize = 15),
                    grid_height = unit(10, "mm"), grid_width = unit(10, "mm"),
                    title_position = "topcenter", 
                    title = "insertion strand", 
                    size = unit(5, "mm"),
                    title_gp = gpar(fontsize = 20, fontface = "bold.italic"), 
                    title_gap = unit(10.5, "mm"),
                    border = "black", background = "white")
lgd_lines = Legend(at = "CIS", 
                   type = "lines", 
                   legend_gp = gpar(col = "black", lwd = 3), 
                   labels_gp = gpar(fontsize = 15, col = "white"),
                   grid_height = unit(10, "mm"), grid_width = unit(10, "mm"),
                   title_position = "topleft",
                   title = "CIS", 
                   size = unit(5, "mm"),
                   title_gp = gpar(fontsize = 20, fontface = "bold.italic"), 
                   title_gap = unit(10.5, "mm"),
                   border = "black", background = "white")
lgd_list_horizontal = packLegend(lgd_links, lgd_points, lgd_lines, direction = "horizontal", gap = unit(30, "mm"))

draw(lgd_list_horizontal, x = unit(.75, "npc"), just = c("center", "center"))
dev.off()