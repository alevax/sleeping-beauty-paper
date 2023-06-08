##
# Pathway Analysis with aREA
# --------------------------

do_pathways_with_aREA <- function( vpmat , threshold = 1.645, out_dir = "")
{
  message(">>> Pathway Analysis with aREA")
  # https://github.com/ctlab/fgsea
  
  # library(devtools)
  # install_github("ctlab/fgsea")
  library(data.table)
  # library(fgsea)
  library(ggplot2)
  
  source("libs/tools/pathway-analysis.R")
  h <- readRDS("data/MSigDB-gene-sets/msigdb-h-as-regulon.rds")
  c2 <- readRDS("data/MSigDB-gene-sets/msigdb-c2-as-regulon.rds")
  # c5 <- readRDS("../vaxtools/data/MSigDB-gene-sets/msigdb-c5-as-regulon.rds")
  # c6 <- readRDS("../vaxtools/data/MSigDB-gene-sets/msigdb-c6-as-regulon.rds")
  # c7 <- readRDS("../vaxtools/data/MSigDB-gene-sets/msigdb-c7-as-regulon.rds")
  
  samples.list <- list(vpmat)
  names(samples.list) <- "vpmat"
  # x <- convertFeaturesOnMatrix(gesmat,"ENTREZID","SYMBOL")
  # samples.list <- list(x)
  # names(samples.list) <- "gesmat"		
  
  reactome_pw <- c2[ grepl("REACTOME",names(c2)) ]
  kegg_pw <- c2[ grepl("KEGG",names(c2)) ]
  
  hormone_pw <- c2[ grepl("ANDRO|ESTRO|GLUCOCO",names(c2)) ]
  
  msigdb.list <- list(h,reactome_pw,kegg_pw)
  names(msigdb.list) <- c("h","reactome","kegg")
  
  system.time({
    pathway.results <- list()
    for( msigdb.index in seq_along(msigdb.list) )
    {
      pathway.results[[names(msigdb.list)[[msigdb.index]]]]
      for( i in seq_along(samples.list) )
      {
        filename <- paste0(out_dir, paste0( names(samples.list)[[i]] , "-" , names(msigdb.list)[[msigdb.index]] ))
        filename.pdf <- paste0(filename,".pdf")
        filename.xls <- paste0(filename,".xls")
        pathway.results[[ names(msigdb.list)[[msigdb.index]] ]][[i]] <- doPathwayAnalysis( samples = samples.list[[i]] , filename = filename.pdf , 
                                                                                           msigdb = NULL , 
                                                                                           threshold = threshold ,
                                                                                           regulon = msigdb.list[[msigdb.index]] ,
                                                                                           save2excel = FALSE , filename.xls = filename.xls , .DEBUG = FALSE )
        # doPathwayAnalysis( samples = samples.list[[i]] , filename = names(samples.list)[[i]] , msigdb = c5 , regulon = genesSetList.as.regulon , threshold = 5 )
      }
    }
  })        
  
  filetag <- "h"
  vpmat_pathways_hallmarks <- unlist(pathway.results,recursive = F)[[filetag]]
  filetag <- "reactome"
  vpmat_pathways_reactome <- unlist(pathway.results,recursive = F)[[filetag]]
  filetag <- "kegg"
  vpmat_pathways_kegg <- unlist(pathway.results,recursive = F)[[filetag]]
  
  return( do.call(rbind,list(vpmat_pathways_hallmarks,vpmat_pathways_reactome,vpmat_pathways_kegg)) )
  
}

do_pathways_with_aREA_with_options <- function( vpmat , threshold = 1.645 , pattern = "*" )
{
  message(">>> Pathway Analysis with aREA")
  library(ggplot2)
  
  source("libs/tools/pathway-analysis.R")
  pw_reg <- readRDS("data/MSigDB-gene-sets/msigdb-msigdb.v7.5.1.symbols.gmt-as-regulon.rds")
  
  samples.list <- list(vpmat)
  names(samples.list) <- "vpmat"
  
  filtered_pw <- pw_reg[ grepl(pattern,names(pw_reg)) ]
  
  msigdb.list <- list(filtered_pw)
  names(msigdb.list) <- c("filtered_pw")
  
  system.time({
    pathway.results <- list()
    for( msigdb.index in seq_along(msigdb.list) )
    {
      pathway.results[[names(msigdb.list)[[msigdb.index]]]]
      for( i in seq_along(samples.list) )
      {
        filename <- paste0( names(samples.list)[[i]] , "-" , names(msigdb.list)[[msigdb.index]] )
        filename.pdf <- paste0(filename,".pdf")
        filename.xls <- paste0(filename,".xls")
        pathway.results[[ names(msigdb.list)[[msigdb.index]] ]][[i]] <- doPathwayAnalysis( samples = samples.list[[i]] , filename = filename.pdf , 
                                                                                           msigdb = NULL , 
                                                                                           threshold = threshold ,
                                                                                           regulon = msigdb.list[[msigdb.index]] ,
                                                                                           save2excel = FALSE , filename.xls = filename.xls , .DEBUG = FALSE )
        # doPathwayAnalysis( samples = samples.list[[i]] , filename = names(samples.list)[[i]] , msigdb = c5 , regulon = genesSetList.as.regulon , threshold = 5 )
      }
    }
  })        
  
  filetag <- "filtered_pw"
  vpmat_pathways_hallmarks <- unlist(pathway.results,recursive = F)[[filetag]]
  
  return( vpmat_pathways_hallmarks )
  
}

