##
# The OncoMatch Function
# ----------------------

# `[` <- function(...) base::`[`(...,drop=FALSE)  

OncoMatch <- function( vpmat_to_test , 
                       vpmat_for_cMRs , 
                       tcm_size = 50 , 
                       both_ways = FALSE ,
                       om_min_threshold = FALSE , 
                       intersect_regs = TRUE ,
                       is_mht_correction = "bonferroni",
                       is_return_best_and_worst_matches = FALSE )
{
  require(crayon)
  require(tidyverse)
  require(viper)
  print_msg_info(">>> Running OncoMatch ...")
  
  if (intersect_regs==TRUE)
  {
    regs_in_common <- intersect(rownames(vpmat_to_test),rownames(vpmat_for_cMRs))
    
    dim(vpmat_for_cMRs)
    vpmat_for_cMRs <- vpmat_for_cMRs[regs_in_common,]
    if("numeric" %in% class(vpmat_for_cMRs)){
      vpmat_for_cMRs <- as.matrix(data.frame(vpmat_for_cMRs))
    }
    dim(vpmat_for_cMRs)
    dim(vpmat_to_test)
    vpmat_to_test <- vpmat_to_test[regs_in_common,]
    if("numeric" %in% class(vpmat_to_test)){
      vpmat_to_test <- as.matrix(data.frame(vpmat_to_test))
    }
    dim(vpmat_to_test)
  }
  
  
  if (both_ways==TRUE)
  {
    om_t <- aREA( vpmat_to_test , generateRegulonObjectFromProteinActivityMatrix(vpmat_for_cMRs,n_top = tcm_size/2 ,justWithOnes = T) , minsize = 0 )$nes 
    om_q <- aREA( vpmat_for_cMRs , generateRegulonObjectFromProteinActivityMatrix(vpmat_to_test,n_top = tcm_size/2 ,justWithOnes = T) , minsize = 0 )$nes 
    
    print_msg_warn("*** Found: " , sum(om_t>30) , " scores above NES=30. Saturating its thresholds ***" )
    om_t <- ifelse(om_t>30,30,om_t)
    om_q <- ifelse(om_q>30,30,om_q)
    
    dim(om_t)
    dim(om_q)
    
    om_q <- t(om_q)
    
    print_msg_warn("*** Found: " , sum( is.na(om_t) ) , " scores as NA. Setting them to ZERO ***")
    print_msg_warn("*** Found: " , sum( is.na(om_q) ) , " scores as NA. Setting them to ZERO ***")
    om_t[ is.na(om_t) ] <- 0
    om_q[ is.na(om_q) ] <- 0
    
    # om_q <- ifelse( om_t > 2 & t(om_q) < 2 , -1*om_q , om_q )
    om <- (om_t+om_q)/2
    
  } else {
    
    # om <- sapply( 1:ncol(vpmat_to_test) , function(i) { print_msg_info(">>> Running sample: " , i ) ; aREA( vpmat_to_test[,i] , generateRegulonObjectFromProteinActivityMatrix(vpmat_for_cMRs,n_top = tcm_size/2 ,justWithOnes = T) , minsize = 1 )$nes  } )
    # rownames(om) <- colnames(vpmat_for_cMRs)
    # colnames(om) <- colnames(vpmat_to_test)
    om <- aREA( vpmat_to_test , generateRegulonObjectFromProteinActivityMatrix(vpmat_for_cMRs,n_top = tcm_size/2 ,justWithOnes = T) , minsize = 1 )$nes
    
    print_msg_warn("*** Found: " , sum( is.na(om) ) , " scores as NA. Setting them to ZERO ***")
    om[ is.na(om) ] <- 0
    print_msg_warn("*** Found: " , sum(om>30) , " scores above NES=30. Saturating its thresholds. max(om): " , max(om) %>% round(1) , "***" )

    om <- ifelse(om>30,30,om)
    # print(table(om>30))
  }
  
  # print(dim(om))
  # print(dim(vpmat_to_test))
  # print(dim(vpmat_for_cMRs))
  # 
  # print(head(rownames(om)))
  # print(head(colnames(om)))
  
  x <- which(om == max(om,na.rm = TRUE), arr.ind = TRUE)
  best_match_test_id <- colnames(om)[ x[,2] ]
  best_match_cMRs_id <- rownames(x)[1]
  flag <- NA
  .my_thres <- 0
  
  while ( is.na(flag) )
  {
    # Idea: for the sample that has the best match, keep iterating slowly increasing the threshold until you find the sample that has the lowest OM score. If can't find it, increase the threhsold until you hit it.
    .my_thres = .my_thres + 0.1
    worst_match_cMRs_id <- om[,best_match_test_id][ om[,best_match_test_id] < .my_thres & om[,best_match_test_id] > -.my_thres ][1]
    flag <- worst_match_cMRs_id
    if (length(worst_match_cMRs_id)==0)
      break ;
  }
  
  print_msg_warn(">>> >> GSEA matching example threshold found: " , .my_thres %>% round(2) )
  
  worst_match_cMRs_id <- names(worst_match_cMRs_id)
  
  # om <- t(om)
  
  if ( is_mht_correction != FALSE )
  {
    print_msg_warn("+++ Correcting for Multiple Hypothesis Testing +++")
    om <- apply( t(om), 2, function(x, adjust) {
      p.adjust(pnorm(x, lower.tail=FALSE), adjust)
      # pnorm(x, lower.tail=FALSE)
    }, adjust=is_mht_correction) 
  } else {
    print_msg_warn("+++ No Multiple Hypothesis Correction +++")
    om <- apply( t(om), 2, function(x, adjust) {
      # p.adjust(pnorm(x, lower.tail=FALSE), adjust)
      pnorm(x, lower.tail=FALSE)
    }, adjust="PIPPO")
  }  
  
  om <- -log10( om )
  
  if ( om_min_threshold != FALSE )
  {
    print_msg_warn("*** Zeroing OM scores below a certain threshold | Not for OL purposes ***")
    om <- ifelse( om < om_min_threshold , 0 , om )  
  }
  
  if ( is_return_best_and_worst_matches == TRUE )
  {
    ret <- list(oncomatch_matrix=om,
                best_match_test_id=best_match_test_id,
                worst_match_cMRs_id=worst_match_cMRs_id,
                best_match_cMRs_id=best_match_cMRs_id)
    return(ret)
  }
  
  # print("hi")
  # print(dim(om))
  # print(head(rownames(om)))
  # print(head(colnames(om)))
  # 
  # print(ncol(vpmat_for_cMRs))
  # print(head(rownames(vpmat_for_cMRs)))
  
  # rownames(om) <- colnames()
  
  return(om)
}

##
# TEST
# ----
# entrez <- sample( entrez.mapping.all$ENTREZID , size = 10 )
## entrez <- entrez.mapping.all$ENTREZID[1:10]
# symbols <- getGeneSymbolFromEntrezId( EntrezIds = entrez )
# entrez.returned <- getEntrezIdFromGeneSymbol( symbols )
# stopifnot( identical( entrez , entrez.returned ) )
#
# match( symbols , entrez2geneSymbol.mapped$SYMBOL )
# match( entrez , entrez2geneSymbol.mapped$ENTREZID )
# which( entrez2geneSymbol.mapped$SYMBOL %in% symbols )
# which( entrez2geneSymbol.mapped$ENTREZID %in% entrez )


generateRegulonObjectFromProteinActivityMatrix <- function( proteinActivityMatrix , n_top = 50 , justWithOnes = TRUE , isSymmetric = TRUE )
{
  generateRegulonElementFromGeneSet <- function(x)
  {
    if( isSymmetric )
    {
      y <- sort(x)
      y <- c(head(y,n_top),tail(y,n_top))
    } else
    {
      y <- sort(abs(x),decreasing = TRUE)
      y <- head(y,2*n_top)
    }
    
    if ( justWithOnes == FALSE )
    {
      tfmode <- y
    }
    else
    {
      tfmode = sign(y)
      # tfmode = ifelse( y < 0 , -1 , y )
      # tfmode = ifelse( y > 0 , +1 , tfmode )
    }
    names(tfmode) <- names(y)
    regulonElement <- list( tfmode = tfmode ,
                            likelihood = rep(1, length(tfmode) )
    )
    return(regulonElement)
  }
  regulonObj <- apply( proteinActivityMatrix , 2 , generateRegulonElementFromGeneSet )
  class(regulonObj) <- "regulon"
  return(regulonObj)
}
print_msg_warn <- function(...) {
  info_msg <- yellow $ bgMagenta $ bold
  text <- list(...)
  # text <- str_c(text,sep = "")
  cat(info_msg(text),"\n",sep = "")
}	
print_msg_info <- function(...) {
  info_msg <- yellow $ bold
  text <- list(...)
  # text <- str_c(text,sep = "")
  cat(info_msg(text),"\n",sep = "")
}		
# Load functions
OncoLoopScoresBinning <- function(mat,color_palette="Oranges", cutoff_max = 30) {
  
  my_max <- round(max(mat),0)
  # set up cut-off values
  # breaks <- c(0,5,10,15,20,25,30,my_max+1)
  breaks <- c(seq(from = 0, to = cutoff_max, by = 5), if(my_max+1 > cutoff_max) {my_max+1} else {NULL}) #c(0,5,10,15,20,25,30,my_max+1)
  
  tags <- c()
  for(i in 1:(length(breaks)-1)){
    tags <- c(tags, paste0("[", breaks[i], "-", breaks[i+1], ")"))
  }
  
  
  # specify interval/bin labels
  # tags <- c("[0-5)","[5-10)", "[10-15)", "[15-20)", "[20-25)", "[25-30)",paste0("[30-",my_max+1,")") )
  # tags <- c("[0-5)","[5-10)", "[10-15)", "[15-20)", "[20-25)", "[25-30)",paste0("[30-",my_max+1,")") )
  tags <- factor(tags)
  
  # x <- RColorBrewer::brewer.pal( length(breaks)-1 , "Set1" )
  x <- RColorBrewer::brewer.pal( length(breaks)-1 , color_palette )
  x[1] <- "#FFFFFF"
  names(x) <- tags
  f1 <- x
  # f1 = colorRamp2( as.integer(tags) , x )
  
  group_tags <- apply(mat, 2 , function(x) {
    res <- cut(x,
               breaks=breaks,
               include.lowest=TRUE,
               right=FALSE,
               labels=tags ,
               # labels=FALSE ,
               ordered_result=TRUE) ;
    # res <- factor(res,levels = tags,ordered = TRUE)
  }) 
  rownames(group_tags) <- rownames(mat)
  # inspect bins
  # summary(group_tags)
  
  return(list(mat=group_tags,colors=f1))
}
