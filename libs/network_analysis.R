#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

suppressPackageStartupMessages(library(igraph))
# suppressPackageStartupMessages(library(shiny))
suppressPackageStartupMessages(library(visNetwork))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(circlize))

source("libs/tools/utils.R")
source("libs/tools/cross_species_utils.R")

# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
# -------------------------------------------------------------------------------------------
# ------------------------------------ LOAD DATA FUNCTION -----------------------------------
# -------------------------------------------------------------------------------------------
# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
getMyTable <- function(filename){
  my_table <- readRDS(filename)
  table(my_table$regulator_type)
  names(my_table)[names(my_table) == "gene_human"] <- "human_gene_name"
  return(my_table)
}

# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
# -------------------------------------------------------------------------------------------
# ----------------------------------- GENERAL HELPER FUNCS ----------------------------------
# -------------------------------------------------------------------------------------------
# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
is.odd <- function(x) x %% 2 != 0
is.even <- function(x) x %% 2 == 0
range_standardize <- function(x, newMin = 0, newMax = 1){
  (x - min(x))/(max(x)-min(x)) * (newMax - newMin) + newMin 
}

# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
# -------------------------------------------------------------------------------------------
# --------------------------------- GET EDGE TABLE FUNCTIONS --------------------------------
# -------------------------------------------------------------------------------------------
# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&

# ------------------------------------- HELPER FUNCTIONS ------------------------------------
# ----------------------------------- GET TOP MR FUNCTIONS ----------------------------------
getTFsWithHighestNEPCScoresDF <- function(my_table, n_top = 50){
  # x = the TFs with the top 50 NEPC VIPER Scores
  top50_TF_df <- my_table %>% 
    filter( regulator_type == "TF" ) %>%
    arrange( NEPC_viper_score ) %>%
    slice_head(n = n_top)
  return(top50_TF_df)
}
getTFsWithLowestNEPCScoresDF <- function(my_table, n_lowest = 50){
  # y = the TFs with the lowest 50 NEPC VIPER Scores
  lowest50_TF_df <- my_table %>% 
    filter( regulator_type == "TF" ) %>%
    arrange( desc(NEPC_viper_score) ) %>%
    slice_head(n = n_lowest)
  return(lowest50_TF_df)
}
getTFsWithTopBottomNEPCScoresDF <- function(my_table, n = 100){
  topN_TF_df <- getTFsWithHighestNEPCScoresDF(my_table, n_top = floor(n/2))
  lowestN_TF_df <- getTFsWithLowestNEPCScoresDF(my_table, n_lowest = floor(n/2))
  # Combine x and y together into z
  topN_lowestN_TF_df <- do.call(rbind,list(topN_TF_df, lowestN_TF_df))
  return(topN_lowestN_TF_df)
}
getTopBottomCandidateMRs <- function(my_table, n = 100){
  topN_lowestN_TF_df <- getTFsWithTopBottomNEPCScoresDF(my_table, n)
  sum( topN_lowestN_TF_df$is_cis_gene == "Yes" | is.na(topN_lowestN_TF_df$NEPC_viper_score) ) # TO FIX
  # Convert the names of the 100 mouse genes into human gene names
  top_mrs <- topN_lowestN_TF_df$human_gene_name# %>% mouse_to_human()
  return(top_mrs)
}

# ------------------------------------- HELPER FUNCTIONS ------------------------------------
# ------------------------------- GET CIS MODULATORS FUNCTIONS ------------------------------
getCandidateModulators <- function(my_table){
  x <- my_table %>% filter(is_cis_gene=="Yes",
                           is_diff_expr=="Yes",
                           is_candidate_modulator=="Yes")
  my_candidate_modulators <- x$human_gene_name
  return(my_candidate_modulators)
}

getCISModulators <- function(my_table){
  my_candidate_modulators <- getCandidateModulators(my_table)
  x <- my_table %>%
    filter( cis_fisher_pvalue < 0.05 ) %>%
    filter( human_gene_name %in% my_candidate_modulators )
  cis_modulators <- x$human_gene_name
  return(cis_modulators)
}
getNEPCMarkers <- function(){
  c(
    "EZH2",
    "ASCL1",
    "ONECUT2",
    "FOXA2",
    # "FOXA1",
    # "SOX2",
    "AR",
    "MYCN",
    "SOX11",
    "FOXM1",
    "CENPF",
    "SYP",
    "CHGA",
    "INSM1",
    "TP53",
    "RB1",
    "PTEN",
    "THBS1",
    # "NKX3-1",
    "KLK3",
    "NEUROD1",
    "POU3F4",
    "SNAI1",
    "TWIST1",
    "TWIST2",
    # "TP53BP1",
    "MYT1"#,
    # "NSD2"
  )
}

# --------------------------------------- MAIN FUNCTION -------------------------------------
getEdgeTable <- function(my_table, cindy_table, n = 100, use_nepc_markers = FALSE){
  if(use_nepc_markers){
    top_mrs <- getNEPCMarkers()
  } else {
    top_mrs <- getTopBottomCandidateMRs(my_table, n)
  }
  cis_modulators <- getCISModulators(my_table)
  edge_table <- cindy_table %>%
    filter(TF %in% top_mrs) %>%
    filter( Modulator %in% cis_modulators)
  edge_table <- edge_table %>% dplyr::rename(from=Modulator,to=TF)
  return(edge_table)
}

# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
# -------------------------------------------------------------------------------------------
# --------------------------------- GET NODE TABLE FUNCTIONS --------------------------------
# -------------------------------------------------------------------------------------------
# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&

# ------------------------------------- HELPER FUNCTIONS ------------------------------------
getNodeTableBasic <- function(edge_table){
  cis_nodes <- unique(edge_table$from)
  mrs_nodes <- unique(edge_table$to)
  node_table <- tibble( id = c(cis_nodes,mrs_nodes) , 
                        type = c(rep("CIS",length(cis_nodes)),rep("MR",length(mrs_nodes)) )
  )
  return(node_table)
}
expandNodeTable <- function(node_table, my_table){
  node_table$gene_name_mouse <- node_table$id %>% human_to_mouse(na.rm = F) 
  node_table$cis_score_size <- round( -log(my_table$cis_fisher_pvalue[ match( node_table$gene_name_mouse , my_table$gene ) ]) / 5 ) + 10
  
  node_table$gene_expr <- my_table$dge_logFC[ match( node_table$gene_name_mouse , my_table$gene ) ]
  node_table$prot_actv <- my_table$NEPC_viper_score[ match( node_table$gene_name_mouse , my_table$gene ) ]
  
  node_table$group <- factor(node_table$type)
  
  return(node_table)
}

# --------------------------------------- MAIN FUNCTION -------------------------------------
getNodeTable <- function(edge_table, my_table){
  node_table <- getNodeTableBasic(edge_table)
  node_table <- expandNodeTable(node_table, my_table)
  return(node_table)
}

# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
# -------------------------------------------------------------------------------------------
# ----------------------------- GET NODE & EDGE TABLE FUNCTIONS -----------------------------
# -------------------------------------------------------------------------------------------
# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
getNodesTableWithAttributesToVisualize <- function(node_table, hide_labels = FALSE){
  # We'll start by adding new node and edge attributes to our dataframes. 
  vis.nodes <- node_table
  # pas_color_fun = colorRamp2(c(-3, 0, 3), c("deepskyblue3", "white", "brown3"))
  pas_color_fun_mrs = colorRamp2(c(-15, -12.5, -10, 0, 10, 12.5, 15),
                             c("deepskyblue4",
                               "deepskyblue",
                               "cadetblue1",
                               "white",
                               "tomato",
                               "red",
                               "red4"))
  pas_color_fun_cis = colorRamp2(#c(-15, -10, 0, 10, 15),
                                 c(-15, 0, 15),
                                 c("deepskyblue",
                                   #"cadetblue1",
                                   "white",
                                   #"lightsalmon",
                                   "tomato2"))
  
  
  #colorRamp2(c(-10, -5, 0, 5, 10), c("deepskyblue4", "deepskyblue3", "white", "brown3"))
  
  vis.nodes$shadow <- TRUE # Nodes will drop shadow
  vis.nodes$title  <- vis.nodes$type # Text on click
  vis.nodes$borderWidth <- 3 # Node border width
  vis.nodes$label  <- vis.nodes$id # Node label
  
  # For MR, use only the pAct (by color scaling)
  vis.nodes$color.background <- pas_color_fun_mrs(node_table$prot_actv)
  vis.nodes$color.background[vis.nodes$type=="CIS"] <- pas_color_fun_cis(node_table$prot_actv[vis.nodes$type=="CIS"])
  vis.nodes$color.border <- pas_color_fun_mrs(node_table$prot_actv)
  # For the CIS, gexpr on the outer and the pAct in the inner (by color scaling)
  vis.nodes$color.border[vis.nodes$type=="CIS"] <- pas_color_fun_cis(node_table$gene_expr[vis.nodes$type=="CIS"])
  # For any CIS or other nodes where the pAct is NA, use the gene expression
  vis.nodes$color.background[is.na(node_table$prot_actv)] <-
    pas_color_fun_cis(node_table$gene_expr[is.na(node_table$prot_actv)])
  
  vis.nodes$color.highlight.background <- "orange"
  vis.nodes$color.highlight.border <- "darkred"
  vis.nodes$font.color <- "black"
  vis.nodes$font.color[vis.nodes$type=="CIS"] <- "black"
  vis.nodes$font.color[vis.nodes$type=="MR"] <- "white"
  vis.nodes$font.background <- NULL
  
  
  # There are two types of nodes. One type has the label inside of it and the other type has the label underneath it. The types with the label inside of it are: ellipse, circle, database, box, text. The ones with the label outside of it are: image, circularImage, diamond, dot, star, triangle, triangleDown, square and icon.
  vis.nodes$shape="circle"
  vis.nodes$shape[vis.nodes$type == "CIS"]="box" #"square"
  vis.nodes$font.size <- 250/nchar(vis.nodes$label)
  vis.nodes$font.size[vis.nodes$type == "CIS"] <- 125 #50
  
  if(hide_labels==TRUE){
    vis.nodes$font.color <- vis.nodes$color.background
  }
  
  # if(hide_labels==TRUE){
  #   # vis.nodes$label[vis.nodes$type=="CIS" & vis.nodes$id!="SIRT1"] <- ""
  #   # vis.nodes$shape[vis.nodes$type=="CIS" & vis.nodes$id!="SIRT1"] <- "square"
  #   # equivalent to: vis.nodes$size[vis.nodes$type=="CIS" & vis.nodes$id!="SIRT1"] <- 75
  #   # since circles and boxes are not affected by size parameter
  #   
  #   # Instead of 
  #   # vis.nodes$shape="dot"
  #   # vis.nodes$shape[vis.nodes$type == "CIS"]="square"
  #   # vis.nodes$size <- 70 #MRs only
  #   # vis.nodes$size[vis.nodes$type == "CIS"]  <- 75
  #   
  #   # Instead of changing the nodes, why not just change the font.color to match the color of the nodes?
  #   
  # } else {
  #   vis.nodes$label  <- vis.nodes$id # Node label
  #   # There are two types of nodes. One type has the label inside of it and the other type has the label underneath it. The types with the label inside of it are: ellipse, circle, database, box, text. The ones with the label outside of it are: image, circularImage, diamond, dot, star, triangle, triangleDown, square and icon.
  #   vis.nodes$shape="circle"
  #   vis.nodes$shape[vis.nodes$type == "CIS"]="box" #"square"
  #   vis.nodes$font.size <- 250/nchar(vis.nodes$label)
  #   vis.nodes$font.size[vis.nodes$type == "CIS"] <- 125 #50
  # }
  
  return(vis.nodes)
}
getEdgeTableWithAttributesToVisualize <- function(edge_table){
  # We'll start by adding new node and edge attributes to our dataframes. 
  vis.links <- edge_table
  vis.links$smooth <- TRUE    # should the edges be curved?
  vis.links$shadow <- FALSE    # edge shadow
  return(vis.links)
}

# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
# -------------------------------------------------------------------------------------------
# ----------------------------- GET MULTI RING LAYOUT FUNCTIONS -----------------------------
# -------------------------------------------------------------------------------------------
# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&

# ------------------------------------- HELPER FUNCTIONS ------------------------------------

get_initial_layout_as_tree <- function(vis.links, vis.nodes){
  graph = graph_from_data_frame(d = vis.links[, c("from", "to")], vis.nodes$id)
  
  layout <- data.frame(layout_as_tree(graph, root=NULL))
  rownames(layout) <- vis.nodes$id
  colnames(layout) <- c("X", "Y")
  layout$node_type <- vis.nodes$type
  layout$protein_activity <- vis.nodes$prot_actv
  layout$gene_expr <- vis.nodes$gene_expr
  layout$cis_score_size <- vis.nodes$cis_score_size
  layout$node_group_by_weight <- vis.nodes$node_group_by_weight
  layout <- layout %>%
    rownames_to_column("id") %>%
    as_tibble() %>%
    arrange(node_type, protein_activity)
  return(layout)
}
round_half_moon <- function(vec) {
  y <- sin(seq(0, pi, length.out = length(vec))) # generate sine values
  round_y <- round(y, 2) * max(vec) # round sine values and scale to match original range of Y coordinates
  return(round_y)
}
round_half_moon_sepPosNegMRs <- function(incr_decr_seq){
  if(length(incr_decr_seq)==1){
    return(round_half_moon(incr_decr_seq))
  }
  if(length(incr_decr_seq)==2){
    incr_seq <- incr_decr_seq[1]
    decr_seq <- incr_decr_seq[2]
  } else {
    midpoint <- which(incr_decr_seq == max(incr_decr_seq))[1]
    incr_seq <- incr_decr_seq[1:midpoint]
    decr_seq <- incr_decr_seq[(midpoint+1):length(incr_decr_seq)]
  }
  
  incr_seq_rounded <- sin(seq(0, pi/2, length.out = length(incr_seq)))
  decr_seq_rounded <- rev(sin(seq(pi, pi/2, length.out = length(decr_seq))))
  
  incr_decr_seq_rounded <- round(c(incr_seq_rounded, decr_seq_rounded)*max(incr_decr_seq),2)
  names(incr_decr_seq_rounded) <- names(incr_decr_seq_rounded)
  stopifnot(length(incr_decr_seq_rounded) == length(incr_decr_seq))
  return(incr_decr_seq_rounded)
}

reposition_CIS_genes_in_straight_line_for_loops <- function(layout){
  layout[is.na(layout$protein_activity), ]$protein_activity <- layout[is.na(layout$protein_activity), ]$gene_expr
  layout <- arrange(layout, protein_activity)
  layout[layout$node_type=="CIS", ]$X <- layout[layout$node_type=="CIS", ]$protein_activity %>%
    rank(na.last = FALSE)
  if(length(layout[layout$node_type=="CIS", ]$X) > 1){
    layout[layout$node_type=="CIS", ]$X <- layout[layout$node_type=="CIS", ]$X %>%
      range_standardize(newMin = -50, newMax = 50)
  }
  # layout[layout$node_type=="CIS", ]$Y <- round_half_moon_sepPosNegMRs(getIncrDecrSeq(n = sum(layout$node_type=="CIS")))*-10 #Archway Pattern
  # layout[layout$node_type=="CIS", ]$Y <- getIncrDecrSeq(n = sum(layout$node_type=="CIS"))*-10 #Archway Pattern
  # layout[layout$node_type=="CIS", ]$Y <- rep(c(1,2),length.out=17)*-10 #Alternating pattern
  layout[layout$node_type=="CIS", ]$Y <- (rep(c(1,2,3),length.out=sum(layout$node_type=="CIS"))*-15) -20 #Alternating pattern
  return(layout)
}
reposition_MRs_in_multi_half_rings <- function(layout){
  node_groups <- sort(na.omit(unique(layout$node_group_by_weight)), decreasing = TRUE)
  n_groups <- length(node_groups)
  n_mrs <- sum(layout$node_type=="MR")
  max_group_size <- max(table(layout$node_group_by_weight))
  y_scaling_factor <- 20
  max_abs_x_coordinate <- 50
  for(i in 1:n_groups){
    group_i <- node_groups[i]
    mr_ids_in_group_i_posAct <- layout %>%
      filter(node_group_by_weight == group_i) %>%
      filter(protein_activity >= 0) %>%
      arrange(protein_activity) %>%
      pull(id)
    mr_ids_in_group_i_negAct <- layout %>%
      filter(node_group_by_weight == group_i) %>%
      filter(protein_activity < 0) %>%
      arrange(protein_activity) %>%
      pull(id)
    mr_ids_in_group_i <- c(mr_ids_in_group_i_posAct, mr_ids_in_group_i_negAct)
    if(length(mr_ids_in_group_i_posAct)!=0){
      layout[layout$id %in% mr_ids_in_group_i_posAct, "X"] <-
        (0 + 1:length(mr_ids_in_group_i_posAct))/max_group_size*max_abs_x_coordinate*2
      layout[layout$id %in% mr_ids_in_group_i_posAct, "Y"] <- i*y_scaling_factor - 1:length(mr_ids_in_group_i_posAct)*2
      # Scale X coordinates
      # Do mathematics to range standardize the X values for each set between set thresholds (-X to +X)
      # Divide up the X range so that the higher group i is the farther out it is positioned
      x_coords_group_i_pos <- layout[layout$id %in% mr_ids_in_group_i_posAct,]$X
      layout[layout$id %in% mr_ids_in_group_i_posAct,]$X <-
        range_standardize(log(abs(x_coords_group_i_pos)+1), #apply log transform to stretch out rings
                          3+i*0.25*(100/n_mrs), #3
                          max_abs_x_coordinate * i/ n_groups)
    }
    if(length(mr_ids_in_group_i_negAct)!=0){
      layout[layout$id %in% mr_ids_in_group_i_negAct, "X"] <-
        (0 - rev(1:(length(mr_ids_in_group_i_negAct))))/max_group_size*max_abs_x_coordinate*2
      layout[layout$id %in% mr_ids_in_group_i_negAct, "Y"] <-
        i*y_scaling_factor - rev(1:length(mr_ids_in_group_i_negAct))*2
      # Scale X coordinates
      x_coords_group_i_neg <- layout[layout$id %in% mr_ids_in_group_i_negAct, ]$X
      if(length(mr_ids_in_group_i_negAct)==1){
        layout[layout$id %in% mr_ids_in_group_i & layout$X < 0,]$X <- 0
      } else {
        layout[layout$id %in% mr_ids_in_group_i_negAct, ]$X <-
          range_standardize(x = log(abs(x_coords_group_i_neg)+1)*-1, #apply log transform to stretch out rings
                            newMin = max_abs_x_coordinate * i/ n_groups*(-1),
                            newMax = -(3+i*0.25)*(100/n_mrs)) #-3
      }
    }
    # # Build a triangle and then round it out
    layout[layout$id %in% mr_ids_in_group_i, "Y"] <-
      round_half_moon_sepPosNegMRs(layout[layout$id %in% mr_ids_in_group_i, ]$Y)+10#+i*10#+5
  }
  return(layout)
}

# -------------------------------------- MAIN FUNCTION -------------------------------------
my_layout_as_multi_half_rings <- function(vis.links, vis.nodes){
  layout <- get_initial_layout_as_tree(vis.links, vis.nodes) %>%
    reposition_MRs_in_multi_half_rings() %>%
    reposition_CIS_genes_in_straight_line_for_loops() %>%
    dplyr::arrange(protein_activity)
  layout_matrix <- layout[match(vis.nodes$id, layout$id),] %>%
    dplyr::select(X,Y) %>%
    as.matrix()
  return(layout_matrix)
}

# -------------------------------------- UNUSED FUNCTIONS -------------------------------------
# separate_MRs_at_top_of_triangle <- function(layout){
#   # Make a little space between the positive and negative MRs
#   layout$X[layout$protein_activity < 0 & layout$node_type == "MR"] <-
#     layout$X[layout$protein_activity < 0 & layout$node_type == "MR"] %>%
#     range_standardize(1, 48.5)
#   layout$X[layout$protein_activity >= 0 & layout$node_type == "MR"] <-
#     layout$X[layout$protein_activity >= 0 & layout$node_type == "MR"] %>%
#     range_standardize(51.5, 100)
#   return(layout)
# }
# getIncrDecrSeq <- function(n){
#   c(1:floor(n/2), if(is.odd(n)){floor(n/2)+1} else {NULL}, floor(n/2):1)
# }
# reposition_MRs_in_triangle <- function(layout){
#   n_mrs <- layout %>%
#     filter(node_type=="MR") %>%
#     pull(id) %>%
#     length()
#   
#   MRs_names_ordered_by_pAct_increasing <- layout %>%
#     filter(node_type=="MR") %>%
#     arrange(protein_activity) %>%
#     pull(id)
#   
#   inc_decr_seq <- getIncrDecrSeq(n_mrs)
#   
#   inc_decr_seq_rounded <- inc_decr_seq
#   inc_decr_seq_rounded <- -(inc_decr_seq_rounded*5)
#   
#   for(i in 1:n_mrs){
#     mr_name <- MRs_names_ordered_by_pAct_increasing[i]
#     layout[layout$id==mr_name, "X"] <- i
#     layout[layout$id==mr_name, "Y"] <- inc_decr_seq_rounded[i]
#   }
#   layout <- separate_MRs_at_top_of_triangle(layout)
#   return(layout)
# }
# reposition_CIS_genes_in_straight_line <- function(layout){
#   layout[layout$node_type=="CIS", ]$X <- layout[layout$node_type=="CIS", ]$protein_activity %>%
#     rank() %>%
#     range_standardize(newMin = 10, newMax = 90)
#   return(layout)
# }
# # Create a custom layout function
# my_layout_as_tree <- function(vis.links, vis.nodes) {
#   # Use layout_as_tree to get the initial layout
#   # graph <- graph_from_data_frame(d = vis.links[, c("from", "to")], vis.nodes$id)
#   # protein_activity = vis.nodes$prot_actv
#   # node_type = vis.nodes$type
#   layout <- get_initial_layout_as_tree(vis.links, vis.nodes) %>%
#     reposition_MRs_in_triangle() %>%
#     reposition_CIS_genes_in_straight_line() %>%
#     dplyr::arrange(protein_activity)
#   layout_matrix <- layout[match(vis.nodes$id, layout$id),] %>%
#     dplyr::select(X,Y) %>%
#     as.matrix()
#   return(layout_matrix)
# }
# my_layout_in_circle <- function(vis.links_sirt1, vis.nodes_sirt1){
#   
#   # vis.nodes_sirt1 <- vis.nodes_sirt1 %>%
#   #   arrange(prot_actv)
#   
#   graph = graph_from_data_frame(d = vis.links_sirt1[, c("from", "to")],
#                                 vis.nodes_sirt1$id)
#   
#   layout <- data.frame(layout_in_circle(graph, order = order(vis.nodes_sirt1$prot_actv)))
#   rownames(layout) <- vis.nodes_sirt1$id
#   colnames(layout) <- c("X", "Y")
#   layout["SIRT1",] <- c(0, 0)
#   
#   return(as.matrix(layout))
# }


# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
# -------------------------------------------------------------------------------------------
# ------------------------------ EDGE AND NODE WEIGHT FUNCTIONS -----------------------------
# -------------------------------------------------------------------------------------------
# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&

# ------------------------------------- HELPER FUNCTIONS ------------------------------------
group_thresholds <- function(data, groupsizes) {
  if (sum(groupsizes) != length(data)) {
    stop(paste0("Group sizes must sum to the length of the data vector (", length(data), ")"))
  }
  data <- sort(data, decreasing = FALSE)
  thresholds <- numeric(length(groupsizes) + 1)
  thresholds[1] <- data[1]
  thresholds[length(thresholds)] <- data[length(data)]
  index <- 0
  for(i in 2:length(groupsizes)){
    index <- index + groupsizes[i-1]
    thresholds[i] <- mean(c(data[index], data[index+1]))
  }
  return(thresholds)
}

# -------------------------------------- MAIN FUNCTIONS -------------------------------------
getEdgeWeight <- function(vis.links){
  # Scale the significant triples for each modulator so we see the top for each group.
  vis.links <- vis.links %>%
    group_by(from) %>%
    mutate(edge_weight = range_standardize(significantTriplets, 1, 100)) %>%
    ungroup()
  return(vis.links)
}
getNodeWeight <- function(vis.nodes, vis.links){
  node_weight_df <- vis.links %>%
    group_by(to) %>%
    summarize("node_weight" = sum(edge_weight))
  vis.nodes <- left_join(vis.nodes, node_weight_df, by = c("id" = "to"))
  vis.nodes$node_weight[vis.nodes$type=="MR"] <- range_standardize(vis.nodes$node_weight[vis.nodes$type=="MR"], 1, 99)
  return(vis.nodes)
}
# Write code so that the nodes are grouped such that there are equal numbers of + and - MRs in the same group.
# Create 2 separate node_weight_group_df and cbind them together
getNodeGroupsByWeight <- function(vis.nodes, pos_group_sizes, neg_group_sizes){
  node_weight_df <- vis.nodes %>%
    filter(type=="MR") %>%
    dplyr::select(id, node_weight, prot_actv)
  node_weight_df_pos <- node_weight_df %>%
    filter(prot_actv >= 0) %>%
    dplyr::select(id, node_weight)
  node_weight_df_neg <- node_weight_df %>%
    filter(prot_actv < 0) %>%
    dplyr::select(id, node_weight)
  pos_thresholds <- group_thresholds(node_weight_df_pos$node_weight, pos_group_sizes)
  neg_thresholds <- group_thresholds(node_weight_df_neg$node_weight, neg_group_sizes)
  node_weight_df_pos$node_group_by_weight <- NA
  node_weight_df_neg$node_group_by_weight <- NA
  for(i in 1:(length(pos_thresholds)-1)){
    node_weight_df_pos$node_group_by_weight[node_weight_df_pos$node_weight >= pos_thresholds[i] & node_weight_df_pos$node_weight < (pos_thresholds[i+1]+0.001)] <- i
  }
  for(i in 1:(length(neg_thresholds)-1)){
    node_weight_df_neg$node_group_by_weight[node_weight_df_neg$node_weight >= neg_thresholds[i] & node_weight_df_neg$node_weight < (neg_thresholds[i+1]+0.001)] <- i
  }
  node_weight_group_df <- rbind(node_weight_df_pos, node_weight_df_neg)
  vis.nodes <- left_join(vis.nodes, node_weight_group_df, by = c("id" = "id"))
  return(vis.nodes)
}

# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
# -------------------------------------------------------------------------------------------
# --------------------------- GET INITIAL NETWORK OBJECT FUNCTION ---------------------------
# -------------------------------------------------------------------------------------------
# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
getNetworkVisualizationObject <- function(vis.nodes, vis.links){
  visnet_1 <- visNetwork(vis.nodes,
                         vis.links,
                         # main = "Integraive Analysis",
                         # submain = "CIS genes as initiation events of NEPC MRs",
                         # width = "500px", height = "400px"
                         #width = "1000px", height = "800px"
                         width = "1500px", height = "1200px"
                         # width = "1550px", height = "1240" #1500*1.2
  ) %>%
    # visGroups(groupname = "CIS",
    #           shape = "box",
    #           color = list(background = "azure", border="black") ) %>%
    # visGroups(groupname = "MR",
    #           shape = "circle",
    #           color = list(background = "gold", border="black") ) %>%
    # Sets global options
    visEdges(arrows = list(to = list(enabled = TRUE, scaleFactor = 0.5, type = 'arrow'))) 
  return(visnet_1)
}

# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
# -------------------------------------------------------------------------------------------
# ------------------------------- VISUALIZATION MAIN FUNCTIONS ------------------------------
# -------------------------------------------------------------------------------------------
# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
getAllCISGenesVisualization <- function(my_table, cindy_table, n_mrs = 100, hide_labels = FALSE){
  require(tidyverse)  	
  require(circlize)  	
  edge_table <- getEdgeTable(my_table, cindy_table, n = n_mrs) %>%
    filter(from != "GPR21") #Filter out GPR21 because its protein activity is NA
  node_table <- getNodeTable(edge_table, my_table)
  
  vis.nodes <- getNodesTableWithAttributesToVisualize(node_table, hide_labels)
  if(n_mrs == 50){
    vis.nodes$font.size[vis.nodes$type == "CIS"] <- 100 #smaller font size when fewer MRs
  }
  
  vis.links <- getEdgeTableWithAttributesToVisualize(edge_table)
  
  # Color the edges and set the width by the # of sig-triplets
  edge_color_fun = colorRamp2(c(0,1), c("#ededed", "darkgray"))
  sirt_edge_color_fun = colorRamp2(c(0,1), c("#FCAE1E", "#ED820E")) # Merigold to Apricot
  
  
  
  vis.links <- getEdgeWeight(vis.links)
  # Distinguish edges from SIRT1 so they can modified to stand out
  vis.links <- mutate(vis.links, is_sirt1 = (vis.links$from=="SIRT1"))
  
  vis.links$color.color=edge_color_fun(range_standardize(vis.links$edge_weight^3))
  vis.links$color.color[vis.links$from=="SIRT1"]=
    sirt_edge_color_fun(range_standardize(vis.links$edge_weight[vis.links$from=="SIRT1"]^3))
  vis.links$width <- 1 #Set width=1 for all links that aren't from SIRT1
  vis.links$width[vis.links$from=="SIRT1"] <-
    range_standardize(log(vis.links$edge_weight[vis.links$from=="SIRT1"], base = 2)^8, 1, 10)
  vis.links <- arrange(vis.links, is_sirt1, edge_weight)
  
  # Scale all edges from 1 to 100
  # For all nodes, give them a "edge-score", which is the sum of the weight of the edges connected to them
  # Using this score (scaled from 1 to 100), separate the nodes into groups (e.g. 1-25,26-50,51-75,76-100)
  # Nodes with the highest score are in the front, lowest score farthest away
  # Keep red nodes on the right and blue on the left
  vis.nodes <- getNodeWeight(vis.nodes, vis.links)
  if(n_mrs == 100){
    vis.nodes <- getNodeGroupsByWeight(vis.nodes,
                                       pos_group_sizes = c(14, 12, 10, 6, 5, 3),
                                       neg_group_sizes = c(14, 12, 10, 6, 5, 3))
  } else if(n_mrs == 50){
    vis.nodes <- getNodeGroupsByWeight(vis.nodes,
                                       pos_group_sizes = c(8, 7, 5, 3, 2),
                                       neg_group_sizes = c(8, 7, 5, 3, 2))
  } else {
    stop("Can only use 50 or 100 MRs. Please set n_mrs=50 or n_mrs=100.")
  }
  
  
  my_layout_as_rings_df <- my_layout_as_multi_half_rings(vis.links, vis.nodes)
  
  # my_igraph = graph_from_data_frame(d = vis.links[, c("from", "to")])#, vis.nodes$id)
  # 
  # cytoscapePing()
  # cytoscapeVersionInfo()
  # my_cytoscape_graph <- RCy3::createNetworkFromIgraph(igraph = my_igraph)
  # my_cytoscape_graph <- RCy3::layoutNetwork(layout.name = "attributes-layout", network = my_cytoscape_graph)#, positions =
  # #my_layout_as_rings_df)
  # class(my_cytoscape_graph)
  # View(my_cytoscape_graph)

  
  
  visnet_1 <- getNetworkVisualizationObject(vis.nodes, vis.links)
  my_visualization <- visnet_1 %>%
    # visLegend() %>%
    visExport(label = "Save as PDF",type = "pdf") %>%
    visIgraphLayout( layout = "layout.norm",
                     randomSeed = 666,
                     smooth = FALSE,
                     physics = FALSE,
                     type = "full",
                     layoutMatrix = my_layout_as_rings_df) %>%
    visPhysics(enabled = FALSE)
  my_visualization
}
getSIRT1NetworkVisualization <- function(my_table,
                                         cindy_table,
                                         n_mrs = 100,
                                         hide_labels = FALSE,
                                         use_nepc_markers = FALSE){
  # cindy_table %>% filter(Modulator=="SIRT1", #TF=="INSM1")
  cindy_table_sirt1 <- filter(cindy_table, Modulator=="SIRT1") %>%
    arrange(desc(significantTriplets))
  
  edge_table <- getEdgeTable(my_table, cindy_table_sirt1, n = n_mrs, use_nepc_markers)# %>%
  node_table <- getNodeTable(edge_table, my_table)
  
  vis.nodes_sirt1 <- getNodesTableWithAttributesToVisualize(node_table, hide_labels)
  vis.links_sirt1 <- getEdgeTableWithAttributesToVisualize(edge_table)
  
  # Color the edges
  # avocado_edge_color_fun = colorRamp2(c(0,0.5,1), c("yellow", "limegreen", "darkgreen"))
  sirt_edge_color_fun = colorRamp2(c(0,1), c("#FFD580", "#ED820E"))
  vis.links_sirt1 <- getEdgeWeight(vis.links_sirt1)
  vis.links_sirt1$color.color=sirt_edge_color_fun(range_standardize(vis.links_sirt1$edge_weight))
  vis.links_sirt1$width=8 #range_standardize(log(vis.links_sirt1$edge_weight, base = 2)^8, 8, 10) #5, 40)
  vis.links_sirt1 <- arrange(vis.links_sirt1, edge_weight)
  
  vis.nodes_sirt1 <- getNodeWeight(vis.nodes_sirt1, vis.links_sirt1)
  
  if(use_nepc_markers==TRUE){
    vis.nodes_sirt1 <- getNodeGroupsByWeight(vis.nodes_sirt1,
                                             pos_group_sizes = c(5, 3, 2), #c(20, 13, 10, 5),
                                             neg_group_sizes = c(3, 2, 1))
    vis.nodes_sirt1$font.size[vis.nodes_sirt1$type=="MR"] <- vis.nodes_sirt1$font.size[vis.nodes_sirt1$type=="MR"]/1.5
    vis.nodes_sirt1$font.color[vis.nodes_sirt1$prot_actv < 8 & vis.nodes_sirt1$prot_actv > -8] <-"#565656" #"#414141" #"#333333"
  } else if(n_mrs == 100){
    vis.nodes_sirt1 <- getNodeGroupsByWeight(vis.nodes_sirt1,
                                             pos_group_sizes = c(14, 11, 9, 6, 5, 3), #c(20, 13, 10, 5),
                                             neg_group_sizes = c(14, 12, 10, 6, 5, 3)) #c(21, 14, 10, 5))
  } else if(n_mrs == 50){
    vis.nodes_sirt1 <- getNodeGroupsByWeight(vis.nodes_sirt1,
                                             pos_group_sizes = c(8, 6, 5, 3, 2), #c(20, 13, 10, 5),
                                             neg_group_sizes = c(9, 6, 5, 3, 2)) #c(21, 14, 10, 5))
  } else {
    stop("Can only use 50 or 100 MRs. Please set n_mrs=50 or n_mrs=100.")
  }
  
  
  my_layout_as_rings_df <- my_layout_as_multi_half_rings(vis.links_sirt1, vis.nodes_sirt1)
  
  
  sirt1_net_visualization <- getNetworkVisualizationObject(vis.nodes_sirt1, vis.links_sirt1) %>%
    # visLegend() %>%
    visExport(label = "Save as PDF", type = "pdf") %>%
    visIgraphLayout( layout = "layout.norm",
                     randomSeed = 666,
                     smooth = FALSE,
                     physics = FALSE,
                     type = "square",
                     layoutMatrix = my_layout_as_rings_df) %>%
    visPhysics(enabled = FALSE)
  sirt1_net_visualization
}

#####################################################################################
#####################################################################################
#################################                   #################################
#################################                   #################################
#################################                   ################################# 
################################# START OF ANALYSIS #################################
#################################                   #################################
#################################                   #################################
#################################                   #################################
#####################################################################################
#####################################################################################

cindy_table <- readRDS("data/cindy-results/cindy-table-su2c-plus-tcga-merged.rds")
my_table <- getMyTable("experiments/integrative-analysis/processed_data/tibble-of-sb-cis-nepc-viper-cindy-integrative-analysis.rds")

reports_dir <- "experiments/network-analysis/reports/"
dir.create(reports_dir, showWarnings = FALSE)
dir.create(paste0(reports_dir, "all_CIS_genes_network_100MRs/"))
dir.create(paste0(reports_dir, "all_CIS_genes_network_50MRs/"))
dir.create(paste0(reports_dir, "sirt1_visualization_100MRs/"))
dir.create(paste0(reports_dir, "sirt1_visualization_50MRs/"))
dir.create(paste0(reports_dir, "sirt1_visualization_NEPCMarkers/"))
# dir.create(paste0(reports_dir, "hidden_labels/"))
# dir.create(paste0(reports_dir, "hidden_labels/sirt1_visualization_50MRs/"))
# dir.create(paste0(reports_dir, "hidden_labels/all_CIS_genes_network_50MRs/"))


# all_CIS_genes_visualization_100MRs <- getAllCISGenesVisualization(my_table, cindy_table, n_mrs = 100)
# saveRDS(all_CIS_genes_visualization_100MRs, "data/all_CIS_genes_visualization_100MRs.rds")

# # Define UI for application that draws a histogram
# ui <- fluidPage(
#   # Application title
#   titlePanel("Integrative Network Analysis"),
#   visNetworkOutput("network_proxy_nodes", width = "500px", height = "375px") #width = "1000px", height = "750px")
# )
# 
# # Define server logic required to draw a histogram
# server <- function(input, output) {
#   
#   output$network_proxy_nodes <- renderVisNetwork({
#     all_CIS_genes_visualization_100MRs
#   })
# }
# 
# # Run the application 
# shinyApp(ui = ui, server = server)


all_CIS_genes_visualization_100MRs <- getAllCISGenesVisualization(my_table, cindy_table, n_mrs = 100)
all_CIS_genes_visualization_50MRs <- getAllCISGenesVisualization(my_table, cindy_table, n_mrs = 50)
sirt1_visualization_100MRs <- getSIRT1NetworkVisualization(my_table, cindy_table, n_mrs = 100)
sirt1_visualization_50MRs <- getSIRT1NetworkVisualization(my_table, cindy_table, n_mrs = 50)
# 
# all_CIS_genes_visualization_50MRs_hiddenLabels <- getAllCISGenesVisualization(my_table,
#                                                                               cindy_table,
#                                                                               n_mrs = 50,
#                                                                               hide_labels = TRUE)
# sirt1_visualization_50MRs_hiddenLabels <- getSIRT1NetworkVisualization(my_table,
#                                                                        cindy_table,
#                                                                        n_mrs = 50,
#                                                                        hide_labels = TRUE)
# 
sirt1_visualization_NEPCMarkers <- getSIRT1NetworkVisualization(
  my_table,
  cindy_table,
  n_mrs = 50,
  hide_labels = FALSE,
  use_nepc_markers = TRUE
)



# saveRDS(all_CIS_genes_visualization_100MRs, "experiments/network-analysis/processed_data/all_CIS_genes_visualization_100MRs_visNetwork.rds")

visSave(all_CIS_genes_visualization_100MRs,
        file = paste0(reports_dir, "all_CIS_genes_network_100MRs/all_CIS_genes_network_100MRs.html"))
visSave(all_CIS_genes_visualization_50MRs,
        file = paste0(reports_dir, "all_CIS_genes_network_50MRs/all_CIS_genes_network_50MRs.html"))

visSave(sirt1_visualization_100MRs,
        file = paste0(reports_dir, "sirt1_visualization_100MRs/sirt1_visualization_100MRs.html"))
visSave(sirt1_visualization_50MRs,
        file = paste0(reports_dir, "sirt1_visualization_50MRs/sirt1_visualization_50MRs.html"))
visSave(sirt1_visualization_NEPCMarkers,
        file = paste0(reports_dir, "sirt1_visualization_NEPCMarkers/sirt1_visualization_NEPCMarkers.html"))



# visSave(all_CIS_genes_visualization_50MRs_hiddenLabels,
#         file = paste0(reports_dir, "hidden_labels/all_CIS_genes_network_50MRs/all_CIS_genes_network_50MRs_hiddenLabels.html"))
# visSave(sirt1_visualization_50MRs_hiddenLabels,
#         file = paste0(reports_dir, "hidden_labels/sirt1_visualization_50MRs/sirt1_visualization_50MRs_hiddenLabels.html"))

# Hide NSD2 label

# # nepc_top_50_candidate_mrs <- readxl::read_xlsx("/Users/AlexanderWang/Desktop/Califano_Lab_Fall_2021/sleeping-beauty-paper/experiments/integrative-analysis/reports/nepc-top-50-candidate-mrs-table.xlsx")
# nepc_top_candidate_mrs_tbl <- readxl::read_xlsx("/Users/AlexanderWang/Desktop/Califano_Lab_Fall_2021/sleeping-beauty-paper/experiments/integrative-analysis/reports/nepc-signatures-table-final-fisher-integration-table.xlsx")
# nepc_top_candidate_mrs_tbl_sort <- nepc_top_candidate_mrs_tbl %>%
#   dplyr::filter(is_c_mrs == TRUE) %>%
#   dplyr::filter(!is.na(gene_human)) %>%
#   dplyr::filter(!is.na(NEPC_viper_score)) %>%
#   dplyr::arrange(desc(NEPC_viper_score)) %>%
#   dplyr::select(gene_human, NEPC_viper_score)
# writexl::write_xlsx(nepc_top_candidate_mrs_tbl_sort , path = "/Users/AlexanderWang/Desktop/nepc_ranked..xlsx")

# Set up a shiny server. Once you have launched the WebApp Via Shiny, then do the SVG pull
