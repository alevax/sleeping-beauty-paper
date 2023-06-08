# LOAD DATA
orthologous.mapping.ensembl_mouse2human <- read.csv("data/translate/mouse2human.csv", row.names = 1) %>%
  as_tibble() %>%
  dplyr::rename(mouse_ens_id = "mouse_ensembl",
                human_ens_id = "human_ensembl")
orthologous.mapping.ensembl_human2mouse <- read.csv("data/translate/human2mouse.csv", row.names = 1) %>%
  as_tibble() %>%
  dplyr::rename(mouse_ens_id = "mouse_ensembl",
                human_ens_id = "human_ensembl")


# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
# -----------------------------------------------------------------------------------
# ---------------------- HUMAN AND MOUSE ORTHOLOGOUS FUNCTIONS ----------------------
# -----------------------------------------------------------------------------------
# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
getHumanOrthologousFromMouseEnsemblIds <- function(x,debug=FALSE) {
  .query <- tibble(mouse_ens_id = x)
  .df <- left_join(.query, orthologous.mapping.ensembl_mouse2human, by="mouse_ens_id")
  .df <- .df %>% 
    group_by(mouse_ens_id) %>%
    dplyr::filter(row_number()==1 )
  
  .df <- left_join(.query,.df,by="mouse_ens_id")
  
  if (!debug)
    res <- .df %>% dplyr::select(human_ens_id) %>% pull()
  else
    res <- .df %>% dplyr::select(mouse_ens_id,human_ens_id)
  
  return(res)
}
getMouseOrthologousFromHumanEnsemblIds <- function(x, debug=FALSE) {
  .query <- tibble(human_ens_id = x)
  .df <- left_join(.query,orthologous.mapping.ensembl_human2mouse,by="human_ens_id")
  .df <- .df %>% 
    group_by(human_ens_id) %>%
    dplyr::filter(row_number()==1 )
  
  .df <- left_join(.query,.df,by="human_ens_id")
  
  if (!debug)
    res <- .df %>% dplyr::select(mouse_ens_id) %>% pull()
  else
    res <- .df %>% dplyr::select(human_ens_id,mouse_ens_id)      
  
  return(res)
}

# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
# -----------------------------------------------------------------------------------
# ---------------------- GENE AND ENSEMBL CONVERSION FUNCTIONS ----------------------
# -----------------------------------------------------------------------------------
# &-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&-&
getGeneSymbolsFromEnsemblId <- function( ids , organism = "human" )
{
  require("AnnotationDbi")
  require("org.Hs.eg.db")
  require("org.Mm.eg.db")
  db = if( organism == "human" ) org.Hs.eg.db else org.Mm.eg.db 
  mapping <- mapIds( db , keys = ids , column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
  message("- Found NA - " , sum(is.na(mapping)) , " total" )
  x <- mapping[ match( ids , names(mapping) ) ]
  x <- as.character(x)
  x <- ifelse( grepl("NULL",x) | is.na(x) , NA , x )
  return( x )
}

getEnsemblIdsFromGeneSymbols <- function( ids , organism = "human" )
{
  require("AnnotationDbi")
  require("org.Hs.eg.db")
  require("org.Mm.eg.db")
  db = if( organism == "human" ) org.Hs.eg.db else org.Mm.eg.db 
  mapping <- mapIds( db , keys = ids , column = "ENSEMBL", keytype = "SYMBOL", multiVals = "first")
  x <- mapping[ match( ids , names(mapping) ) ]
  x <- as.character(x)
  x <- ifelse( grepl("NULL",x) | is.na(x) , NA , x )
  return( x )
}	



#####################################################################################
#####################################################################################
#################################                   #################################
#################################                   #################################
#################################                   ################################# 
#################################   MAIN FUNCTION   #################################
#################################                   #################################
#################################                   #################################
#################################                   #################################
#####################################################################################
#####################################################################################
human_to_mouse <- function( mat , isEns = FALSE , na.rm = TRUE )
{
  if (is.matrix(mat))
  {
    if ( isEns ) {
      rownames(mat) <- getHumanOrthologousFromMouseEnsemblIds( rownames(mat) )
    } else {
      rownames(mat) <- getGeneSymbolsFromEnsemblId( getMouseOrthologousFromHumanEnsemblIds( getEnsemblIdsFromGeneSymbols( rownames(mat) , "human" ) ) , "mouse")
    }
    
    if ( na.rm )
    {
      res <- mat[ !is.na(rownames(mat)) , ]
      message(">>> Lost " , nrow(mat)-nrow(res) , " features during conversion (NAs)" )        
    }
    
    return(res)
  } else if ( is.vector(mat) & !is.null(names(mat)) ) {
    
    if ( isEns ) {
      names(mat) <- getMouseOrthologousFromHumanEnsemblIds( names(mat) )
    } else {
      names(mat) <- getGeneSymbolsFromEnsemblId( getMouseOrthologousFromHumanEnsemblIds( getEnsemblIdsFromGeneSymbols( names(mat) , "human" ) , "mouse") )
    }
    if ( na.rm )
    {
      res <- mat[ !is.na(names(mat)) ]
      message(">>> Lost " , length(mat)-length(res) , " features during conversion (NAs)" )
    }
    return(res)      
    
  } else if ( is.vector(mat) & is.null(names(mat)) ) {
    
    if ( isEns ) {
      res <- getMouseOrthologousFromHumanEnsemblIds( mat )
    } else {
      res <- getGeneSymbolsFromEnsemblId( getMouseOrthologousFromHumanEnsemblIds( getEnsemblIdsFromGeneSymbols( mat , "human" ) ) , "mouse")        
    }
    
    if ( na.rm )
    {
      res <- res[ !is.na(res) ]
      message(">>> Lost " , length(mat)-length(res) , " features during conversion (NAs)" )
    }
    return(res)        
  } else {
    message("*** ERROR | You must pass vector or matrix as input ***" )
    return( NULL )
  }
  
}  
mouse_to_human <- function( mat , isEns = FALSE , na.rm = TRUE )
{
  if (is.matrix(mat))
  {
    if ( isEns ) {
      rownames(mat) <- getHumanOrthologousFromMouseEnsemblIds( rownames(mat) )
    } else {
      rownames(mat) <- getGeneSymbolsFromEnsemblId( getHumanOrthologousFromMouseEnsemblIds( getEnsemblIdsFromGeneSymbols( rownames(mat) , "mouse" ) ) )
    }
    
    if (na.rm)
    {
      res <- mat[ !is.na(rownames(mat)) , ]
      message(">>> Lost " , nrow(mat)-nrow(res) , " features during conversion (NAs)" )
      message(">>> >> Removing NAs ..." )
    } else {
      res <- mat
      message(">>> >> Keeping NAs ..." )
    }
    return(res)
    
  } else if ( is.vector(mat) & !is.null(names(mat)) ) {
    
    if ( isEns ) {
      names(mat) <- getHumanOrthologousFromMouseEnsemblIds( names(mat) )
    } else {
      names(mat) <- getGeneSymbolsFromEnsemblId( getHumanOrthologousFromMouseEnsemblIds( getEnsemblIdsFromGeneSymbols( names(mat) , "mouse" ) ) )
    }
    
    if (na.rm)
    {
      res <- mat[ !is.na(names(mat)) ]
      message(">>> Lost " , length(mat)-length(res) , " features during conversion (NAs)" )        
      message(">>> >> Removing NAs ..." )
    } else {
      res <- mat
      message(">>> >> Keeping NAs ..." )
    }
    
    return(res)      
    
  } else if ( is.vector(mat) & is.null(names(mat)) ) {
    
    if ( isEns ) {
      res <- getHumanOrthologousFromMouseEnsemblIds( mat )
    } else {
      res <- getGeneSymbolsFromEnsemblId( getHumanOrthologousFromMouseEnsemblIds( getEnsemblIdsFromGeneSymbols( mat , "mouse" ) ) )        
    }
    
    if (na.rm)
    {
      message(">>> >> Removing NAs ..." )
      tmp <- res
      res <- res[ !is.na(res) ]
      message(">>> Lost " , length(tmp)-length(res) , " features during conversion (NAs)" )        
    } else {
      message(">>> >> Keeping NAs ..." )
      res <- res
    }      
    
    return(res)        
    
  } else {
    message("*** ERROR | You must pass vector or matrix as input ***" )
    return( NULL )
  }
  
}




