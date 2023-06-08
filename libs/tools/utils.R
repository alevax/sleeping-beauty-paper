print_msg_info <- function(...) {
  require(crayon)
  info_msg <- yellow $ bold
  text <- list(...)
  # text <- str_c(text,sep = "")
  cat(info_msg(text),"\n",sep = "")
}
print_msg_warn <- function(...) {
  require(crayon)
  info_msg <- black $ bgYellow $ bold
  text <- list(...)
  # text <- str_c(text,sep = "")
  cat(info_msg(text),"\n",sep = "")
}	
init_pheatmap_color_palette <- function( z.limit = 64 , z.satur = 8 , z.signific = 1.6449 , z.lengthout = 7 , color = "RdBu" )
{
  palette.breaks <- c( -z.limit , seq(-z.satur,-z.signific,length.out = z.lengthout ) , 0 , seq(z.signific,z.satur,length.out = z.lengthout ) , +z.limit )
  # color.palette <- colorRampPalette( rev(brewer.pal(11, "RdBu")) )(length(palette.breaks)-1)
  color.palette <- colorRampPalette( rev(brewer.pal(11, color )) )(length(palette.breaks))
  color.palette[c(8,9)] <- "#FFFFFF"
  length(palette.breaks)
  length(color.palette)
  ret <- list(palette.breaks,color.palette)
  names(ret) <- c("palette.breaks","color.palette")
  return( ret )
}
doStouffer <- function( aMatrix , weights = NULL ) {
  if ( is.null(weights) )
  {
    res <- rowSums(aMatrix) / sqrt( ncol(aMatrix) )
  }
  else {
    
    res <- (aMatrix %*% weights) / sqrt( sum(weights**2) )
    # res <- (aMatrix %*% weights) / sqrt( ncol(aMatrix) )
    # res2 <- rowSums(aMatrix * weights) / sqrt( ncol(aMatrix) )
  }
  return(res)
}	
doStoufferFromClusters <- function( mat , clusters , weights = NULL ) 
{
  # colnames(vpmat)
  # identical( colnames(vpmat) , clustering.tibble$sample_id[ match( colnames(vpmat) , clustering.tibble$sample_id ) ] )
  # clusters <- clustering.tibble$cluster_id[ match( colnames(vpmat) , clustering.tibble$sample_id ) ]
  # weights <- clustering.tibble$sil_score[ match( colnames(vpmat) , clustering.tibble$sample_id ) ]
  # mat <- vpmat
  
  if ( !is.null(weights) ) ## If Stouffer's has to be computed without weights
  {
    ## weights : has to be a named vector
    if ( is.null(names(weights)) )
    {
      stop("*** WARNING | must supply weights as named vector")
      return(NA)
    } else {
      stopifnot( identical( colnames(mat) , names(weights) ) )
      stopifnot( identical( colnames(mat) , names(clusters) ) )
    }
    isWithWeights <- TRUE
    print_msg_info(">>> >> Applying Stouffer's integrations with WEIGHTS")    
  } else {
    isWithWeights <- FALSE
    print_msg_info(">>> >> Applying Stouffer's integrations with NO weights")
  }
  
  my_list <- list()
  for( a_clust_id in seq_len(nlevels(clusters)) )
  {
    index <- colnames(mat) %in% names(clusters[clusters == levels(clusters)[a_clust_id]])
    print(sum(index))
    
    if (isWithWeights)
      my_list[[a_clust_id]] <- doStouffer( aMatrix = mat[ , index ] ,
                                           weights = weights[ index ] )
    else
      my_list[[a_clust_id]] <- doStouffer( aMatrix = mat[ , index ] )
  }
  # str(my_list,1)
  names(my_list) <- seq_len(nlevels(clusters))
  mat_stouffered <- do.call(cbind,my_list)	
  colnames(mat_stouffered) <- levels(clusters)
  return( mat_stouffered )
}
rankNorm <- function(x,FUN=median,trim=0)
{
  .rank <- apply(x, 2, rank)
  .median <- apply(.rank, 1, FUN, trim = trim )
  .mad <- apply(.rank, 1, mad)
  x <- (.rank - .median)/.mad
  
  message("- Number of NA features: " , sum( rowSums(is.na(x)) ) ) 
  message("- Number of Inf features: " , sum( rowSums(x) == Inf ) ) 
  message("- Number of 0 features: " , sum( rowSums(x) == 0 ) ) 
  
  message("- Features to Remove:")
  .null.features <- x[ which( rowSums( apply(x,2,is.na) ) > 0 ) , ]
  # print(.null.features)
  message("- Removing NULL/NA features ...")
  x <- na.omit(x)
  message("- Number of NA features: " , sum( rowSums(is.na(x)) ) )
  message("- Number of Inf features: " , sum( rowSums(x) == Inf ) )
  message("- Number of 0 features: " , sum( rowSums(x) == 0 ) ) 
  
  return(x)
}  