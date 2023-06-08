  ##
  # Pathway Analysis - Library
  # --------------------------
	library(pheatmap) # to plot heatmap
  library(viper) # For aREA function

		##
		# GENE SET ANALYSIS
		# -----------------

	  # message(">>> Loading gene sets from MSigDB ...") # Downloaded from http://bioinf.wehi.edu.au/software/MSigDB/index.html
	  #   load( file = "../vaxtools/data/MSigDB-gene-sets/human_c1_v5p2.rdata" )
	  #   load( file = "../vaxtools/data/MSigDB-gene-sets/human_c2_v5p2.rdata" )
	  #   load( file = "../vaxtools/data/MSigDB-gene-sets/human_c3_v5p2.rdata" )
	  #   load( file = "../vaxtools/data/MSigDB-gene-sets/human_c4_v5p2.rdata" )
	  #   load( file = "../vaxtools/data/MSigDB-gene-sets/human_c5_v5p2.rdata" )
	  #   load( file = "../vaxtools/data/MSigDB-gene-sets/human_c6_v5p2.rdata" )
	  #   load( file = "../vaxtools/data/MSigDB-gene-sets/human_c7_v5p2.rdata" )
	  #   load( file = "../vaxtools/data/MSigDB-gene-sets/human_H_v5p2.rdata" )
    #
	  # geneSetCollectionObj <- vector("list", 8 )
	  # geneSetCollectionObj[[1]] <- as.list(Hs.c1)
	  # geneSetCollectionObj[[2]] <- as.list(Hs.c2)
	  # geneSetCollectionObj[[3]] <- as.list(Hs.c3)
	  # geneSetCollectionObj[[4]] <- as.list(Hs.c4)
	  # geneSetCollectionObj[[5]] <- as.list(Hs.c5)
	  # geneSetCollectionObj[[6]] <- as.list(Hs.c6)
	  # geneSetCollectionObj[[7]] <- as.list(Hs.c7)
	  # geneSetCollectionObj[[8]] <- as.list(Hs.H)
    #
		# names(geneSetCollectionObj) <- c( "C1 Positional Gene Sets" , "C2 curated gene sets" , "C3 motif gene sets" , "C4 computational gene sets" ,
		# 																	"C5 GO gene sets" , "C6 oncogenic signatures" , "C7 immunologic signatures" , "H Hallmark gene sets" )
    #
		# str(geneSetCollectionObj,1)

				# x <- read.csv2( "sources/shiny-app/data/hallmarks_Drake_Paull.csv" , stringsAsFactors = FALSE ,header = F ,sep = "\t" )
				# x <- data.frame( t(x) , stringsAsFactors = FALSE )
				# colnames(x) <- x[1,]
				# x <- x[-1,]
				#
				# x <- lapply( x , function(x) x[ x != "" ] )
				# x <- lapply( x , getEntrezIdFromGeneSymbol )
				# saveRDS( object = x , file = "../vaxtools/data/hallmarks-cancer-paull-drake.rds")

	# geneSetCollectionObj[["Hallmark - Paull & Drake 2016"]] <- readRDS( file = "../vaxtools/data/hallmarks-cancer-paull-drake.rds" )
	# str(geneSetCollectionObj,1)

	##
	# Function that perform the enrichment with aREA
	# ----------------------------------------------
		doEnrichment <- function( signature , geneSet , isWithNES = TRUE )
	  {
		  genesSetList.as.regulon <- generateRegulonObjectFromGeneSetList( geneSetList = geneSet )
		  if (isWithNES)
			  tmp <- aREA( signature , genesSetList.as.regulon , minsize = 10 )$nes
		  else
		    tmp <- aREA( signature , genesSetList.as.regulon , minsize = 10 )$es
		  
	    res <- as.vector( tmp )
	    names(res) <- rownames( tmp )
	    return(res)
		}

	##
	# Get the gene set selected by the argument
	# -----------------------------------------
		getGeneSets <- function( setNames = c("hallm") )
		{
			if ( length(setNames) > 1 )
			{
				setNames <- paste0( setNames , collapse = "|" )
			}
			res <- geneSetCollectionObj[ grep( pattern = setNames , x = names(geneSetCollectionObj) , ignore.case = TRUE ) ]
			return(res)
		}


	##
	# Test and examples
	# -----------------

		# 	enrichment.list <- list()
		#   .vpmat <- organoids$pas.expmat.with.tcga
		#   rownames(.vpmat) <- getEntrezIdFromGeneSymbol(rownames(.vpmat))
		#
		#   geneSet.list <- getGeneSets("hallm")
		#
		#   for( index in 1:ncol(.vpmat) )
		#   {
		#     signature <- .vpmat[,index]
		#     sample_name <- colnames(.vpmat)[index]
		#     enrichment.list[[index]] <- doEnrichment( signature = signature , geneSet = geneSet.list[[1]] )
		#   	names(enrichment.list)[index] <- sample_name
		#   }
		#   str(enrichment.list,1)

    ##
    # Pathway Analysis with aREA and MSigDB
    # -------------------------------------
			# require(clusterProfiler)
			generateRegulonsForMSigDB <-  function( gmt.dir = "../vaxtools/data/MSigDB-gene-sets/" , output.dir = "../vaxtools/data/MSigDB-gene-sets/" )
			{
				gmt.file.list <- list.files( path = gmt.dir , pattern = "*.symbols.gmt$" , full.names = FALSE )
				for ( aGmtFile in gmt.file.list )
				{
					geneSetName <- sub( x = aGmtFile , pattern = "(.*).all(\\..*)$" , replacement = "\\1")
					message( "- Generating regulon object for " , geneSetName , " Gene Set ...")
				
					gmtfile <- file.path( gmt.dir , aGmtFile )
					# msigdb <- read.gmt(gmtfile)
					msigdb <- GSA::GSA.read.gmt(gmtfile)
	        .geneSetList <- lapply( levels(msigdb$ont) , function(i) { tmp <- msigdb[msigdb$ont %in% i,]$gene } )
	        names(.geneSetList) <- levels(msigdb$ont)
	        genesSetList.as.regulon <- generateRegulonObjectFromGeneSetList( geneSetList = .geneSetList )				
	        saveRDS( object = genesSetList.as.regulon , file = file.path( output.dir , paste0("msigdb-", geneSetName , "-as-regulon.rds" ) ) )
				}
				message( "- [complete]")
			}
			# system.time(generateRegulonsForMSigDB())
	
      doPathwayAnalysis <- function( samples , msigdb = NULL , threshold = 1.644 , regulon = NULL , filename , .DEBUG = TRUE , save2excel = FALSE , filename.xls = NULL , plot = TRUE ,verbose = FALSE )
      {
        message( ">>> Pathway Analysis with aREA and MSigDB ...")

	        if (is.null(regulon))
	        {
	          message( "- Please wait the building of a regulon object for the analysis ...")
	          .geneSetList <- lapply( levels(msigdb$ont) , function(i) { tmp <- msigdb[msigdb$ont %in% i,]$gene } )
	          names(.geneSetList) <- levels(msigdb$ont)
	          genesSetList.as.regulon <- generateRegulonObjectFromGeneSetList( geneSetList = .geneSetList )
	          
	          enriched.matrix <- matrix( 1 , nrow = length(levels(msigdb$ont)) , ncol(samples) , dimnames = list(levels(msigdb$ont),colnames(samples)) )
	          
	        } else {
	          message( "- Regulon Object already build - skipping ...")
	          genesSetList.as.regulon <- regulon
	          
	          enriched.matrix <- matrix( 1 , nrow = length(genesSetList.as.regulon) , ncol(samples) , dimnames = list(names(genesSetList.as.regulon),colnames(samples)) )
	        }

        message( ">>> Pathway Analysis" )
        {
          pb <- txtProgressBar(min = 0, max = ncol(samples), style = 3)
          for( index in seq_len(ncol(samples)) )
          {
            x <- samples[,index]
            .drug_name <- colnames(samples)[index]
            
            if(verbose) message(">>> Pathway Analysis of " , .drug_name )
            
            res <- aREA( x , genesSetList.as.regulon , minsize = 2 )$nes
            # res <- aREA( x , genesSetList.as.regulon , minsize = 2 )$es
            tmp <- as.vector(res)
            names(tmp) <- rownames(res)
       
            enriched.matrix[ names(tmp) , .drug_name] <- tmp
            
            setTxtProgressBar(pb, index)
          }
          close(pb)
          
          signif.keys.at.least.one.sample <- !(colSums( apply(enriched.matrix,1, function(x) abs(x)<threshold)) == ncol(enriched.matrix) )

          if(.DEBUG)
          {
            a <- apply( enriched.matrix[!signif.keys.at.least.one.sample,] , 1 , min )
            b <- apply( enriched.matrix[!signif.keys.at.least.one.sample,] , 1 , max )
            if( length( a[ a > threshold ] ) != 0 || length( b[ b > threshold ] ) != 0 )
              return("ERROR")
            print(a)
            print(b)
          }

          enriched.matrix <- enriched.matrix[ signif.keys.at.least.one.sample , ]

          if (plot == TRUE)
          {
	          graphics <- init_pheatmap_color_palette( z.limit = 32 , z.signific = threshold )
	          pheatmap( file = filename ,
	                    main = paste0( "Showing " , nrow(enriched.matrix) , "/" , length(genesSetList.as.regulon) , " gene sets" ) ,
	                    mat = enriched.matrix ,
	                    # color = colorRampPalette(rev(brewer.pal(n = 8, name = "BuGn")))(100) ,
	            	      color = graphics$color.palette ,
	                	  breaks = graphics$palette.breaks ,
	          					# clustering_distance_rows = "correlation",
											# clustering_distance_cols = "correlation",
	          					clustering_method = "ward.D2" ,
	                    display_numbers = TRUE ,
	                    number_format = "%.1f" ,
	                    # number_color = "gray" ,
	        	          number_color = "white" ,
	          					# border_color = "gray" ,
	          	border_color = "white" ,
	                    cellwidth = 10 , cellheight = 10 ,
	        						fontsize = 5)
          }
        }

        if ( save2excel )
        {
          library(xlsx)
          x <- as.data.frame( round( enriched.matrix , digits = 3 ) )
          write.xlsx( x = x , file = filename.xls ,
            # sheetName = "Protein Activity Signatures - 6 Hours" ,
            row.names = TRUE , col.names = TRUE )
        }

        return(enriched.matrix)
      }

      