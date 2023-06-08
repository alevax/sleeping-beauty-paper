# source("sources/utils/init.R")
generateRegulonElementFromGeneSet <- function(x)
{
  if ( is.null(names(x)) | class(x) == "character" ) # Create a normal gene set as in the classical GSEA
  {
    tfmode <- as.numeric(rep(1, length(x) ))
    names(tfmode) <- x
  }
  else # Create a special gene set as for aREA
  {
    tfmode <- x
    names(tfmode) <- names(x)
    
    if ( rm.NA == TRUE )
    {
      na.found <- !is.na( names(tfmode) )
      message( " *** Removing " , sum(na.found) , " <NA> elements from the gene set ***")
      tfmode <- tfmode[ na.found ]
    }
  }
  
  regulonElement <- list( tfmode = tfmode ,
                          likelihood = rep(1, length(tfmode) )
  )
  return(regulonElement)
}
generateRegulonObjectFromGeneSetList <- function( geneSetList , rm.NA = FALSE )
{
  regulonObj <- lapply( geneSetList , generateRegulonElementFromGeneSet )
  return(regulonObj)
  # regulonObj <- generateRegulonElementFromGeneSet( geneSetList )
  # names(regulonObj) <-  "The TF"
  # interactomeObj <- list(regulonObj)
  # class(interactomeObj) <- "regulon"
  # return(interactomeObj)
}

message("- Analysis of the NET signature ")
{
  net.signature.table <- read.table( 
    # file.path( "../prostate-cancer/data/from-prem/net-signature.csv" ) 
    file.path( "data/signatures/net-signature.csv" ) 
    , dec = "." , sep = "," , header = T , stringsAsFactors = F )
  # str(net.signature.table)
  net.signature.mrs.all <- net.signature.table$GeneID
  net.signature.mrs.validated <- net.signature.table$Symbol[ net.signature.table$validated == 1 ]
  
  # net.signature.mrs.validated.as.list <- list(getGeneSymbolsFromEntrezId(net.signature.mrs.validated))
  net.signature.mrs.validated.as.list <- list(net.signature.mrs.validated)
  names(net.signature.mrs.validated.as.list) <- "net-signature"
  net.regulon.obj <- generateRegulonObjectFromGeneSetList( net.signature.mrs.validated.as.list )
  class(net.regulon.obj) <- "regulon"
  
  ## Beltran NEPC
  beltran.net.signature.table <- read.table( 
    # file.path( "../prostate-cancer/data/datasets/beltran/beltran-nepc-signature.csv" ) ,
    file.path( "data/signatures/beltran-nepc-signature.csv" ) ,
    dec = "." , sep = "," , header = T , stringsAsFactors = F , as.is = T , skip = 1 )
  # str(beltran.net.signature.table)
  a <- beltran.net.signature.table$HGNC.ID[ beltran.net.signature.table$RNA..CRPC.NE.vs.CRPC.Adeno. == "Over-expressed" ]
  b <- beltran.net.signature.table$HGNC.ID[ beltran.net.signature.table$RNA..CRPC.NE.vs.CRPC.Adeno. == "Under-expressed" ]
  tmp <- c( rep(1,length(a)) , rep(-1,length(b)) )
  names(tmp) <- c(a,b)
  regulonElement <- list( tfmode = tmp ,
                          likelihood = rep(1, length(tmp))
  )
  beltran.net.regulon.obj <- list( regulonElement )
  class(beltran.net.regulon.obj) <- "regulon"
  
  beltran.nepc.and.prad.table <- read.csv2( 
    file.path( "data/signatures/beltran-2016-table-7.csv" ), 
    dec = "." , sep = "," , header = T , stringsAsFactors = F , as.is = T , skip = 1 )
  beltran.nepc.vs.adeno.ges <- beltran.nepc.and.prad.table$CRPC.NE.mean.expression - beltran.nepc.and.prad.table$CRPC.Adeno.mean.expression
  names(beltran.nepc.vs.adeno.ges) <- beltran.nepc.and.prad.table$Hugo_Symbol
  
  beltran.nepc.vs.adeno.ges <- qnorm( beltran.nepc.and.prad.table$P.val.expression/2 , lower.tail = FALSE ) * sign(beltran.nepc.vs.adeno.ges)
  
  # beltran.nepc.vs.adeno.ges[ names(beltran.nepc.vs.adeno.ges) %in% "WWTR1" ]
  
  # beltran.nepc.and.prad.table$P.val.expression
  # beltran.nepc.and.prad.table$P.val.expression[ names(beltran.nepc.vs.adeno.ges) %in% "WWTR1" ]
  
}

getNETenrichment <- function( dataset , regulonObject = net.regulon.obj , score = "nes" )
{
  .score <- match.arg( score , c("nes","es") )
  if ( .score == "es" )
    x <- apply( dataset , 2 , function(x) aREA( x , regulonObject , minsize = 5 )$es )
  else if ( .score == "nes" )
    x <- apply( dataset , 2 , function(x) aREA( x , regulonObject , minsize = 5 )$nes )
  
  return( x )
}

message("- Hieronymus, et al. Cancer Cell (2008) [27 genes]")
{
  # hieronymus_AR_genes <- c("PSA","TMPRSS2","NKX3-1","KLK2","GNMT","TMEPAI","MPHOS9","ZBZB10","EAF2","BM039",
  # 												 "SARG","ACSL3","PTGER4","ABCC4","NNMT","ADAM7","FKBP5","ELL2","MED28","HERC3","MAF",
  # 												 "TNK1","GLRA3","MAPRE2","PIP5K2B","MAN1A1","CD200")
  # hieronymus_AR_regulon <- generateRegulonElementFromGeneSet(hieronymus_AR_genes)
  hieronymus_AR_genes_pos <- c("PSA","KLK3","TMPRSS2","NKX3-1","NKX3.1","KLK2","GNMT","TMEPAI","MPHOS9","ZBZB10","EAF2","BM039",
                               "SARG","ACSL3","PTGER4","ABCC4","NNMT","ADAM7","FKBP5","ELL2","MED28","HERC3","MAF")
  hieronymus_AR_genes_neg <- c("TNK1","GLRA3","MAPRE2","PIP5K2B","MAN1A1","CD200")
  tmp <- c( rep(1,length(hieronymus_AR_genes_pos)) , rep(-1,length(hieronymus_AR_genes_neg)) )
  names(tmp) <- c(hieronymus_AR_genes_pos,hieronymus_AR_genes_neg)
  ar.regulon.obj <- list( hieronymus_AR_genes = list( tfmode = tmp ,
                                                      likelihood = rep(1, length(tmp)))
  )
  class(ar.regulon.obj) <- "regulon"			
  
  # x <- list(hieronymus_AR_genes)
  # names(x) <- "ar-signature"
  # ar.regulon.obj <- generateRegulonObjectFromGeneSetList( x )
  # class(ar.regulon.obj) <- "regulon"			
}
