differentialExpression <- function(dataSetSeries, fileWithGenes) {
  library(GEOquery)
  library(limma)
  library(org.Mm.eg.db)
  library(org.Rn.eg.db)
  library(org.Hs.eg.db)
  library(gtools)
  library(stringr)
  library(fgsea)
  library(ggplot2)
  library(data.table)
  # for collapseBy
  source("https://raw.githubusercontent.com/assaron/r-utils/master/R/exprs.R")
  # for conditions handling
  source("./dif_expression/couple_universal_handler.R")
  # for gene set enrichment analysis
  source("./dif_expression/gsea.R")
  
  # options(download.file.method.GEOquery = "libcurl")
  gse <- getGEO(dataSetSeries, destdir = "./data")[[1]] 

  characteristics <- getCharacteristicsColumns(gse)
  conditionLists <- list()
  explanatoryTable <- data.table()
  
  if (length(characteristics) > 0) {
    message("There are characteristics columns in the samples table.")
    conStructure <- getConditionsFromCharacteristics(gse, characteristics)
    conditionLists <- conStructure$conditionsList
    explanatoryTable <- conStructure$explanatoryTable 
  }
  
  if (length(conditionLists) > 0) {
    pData(gse)$condition <- fillGseConditionColumn(conditionLists)
    message("Column 'condition' was created and filled in the samples table.")
    
    a_gse <<- gse
    resValidation <- provideValidOfSomeColumns(gse)
    fData(gse) <- resValidation$fData
    isContainedIdGenesColumn <- resValidation$successValidation
    
    if (isContainedIdGenesColumn) {
      es <- ""
      if (length(fData(gse)$ENTREZ_GENE_ID) != 0)
        es <- collapseBy(gse, fData(gse)$ENTREZ_GENE_ID, FUN=median)
      else if (length(fData(gse)$GENE) != 0)
        es <- collapseBy(gse, fData(gse)$GENE, FUN=median)

      es <- es[!grepl("///", rownames(es)), ]
      es <- es[rownames(es) != "", ]
      
      message("Garbage was deleted from gene table.")
      a_es00 <<- es
      fData(es) <- data.frame(row.names = rownames(es))
      a_es1 <<- es
      database <- getDatabaseForMapping(gse)
      fData(es)$symbol <- mapIds(database, keys = rownames(es), column = "SYMBOL", keytype = "ENTREZID")
      a_es <<- es
      exprs(es) <- normalizeBetweenArrays(log2(exprs(es)+1), method="quantile")
      
      es.design <- model.matrix(~0+condition, data=pData(es))
      
      fit <- lmFit(es, es.design)
      conditions <- getConditionsForBuildingLinearModel(pData(gse)$condition)
      a_conditions <<- conditions
      
      message("Conditions combinations were received for filling of contrast matrix.")
      
      # Show received pairs of comparisons
      for (i in 1:length(conditions)) 
        message(conditions[[i]][1], " ", conditions[[i]][2])
      
      deSize <- dim(fData(es))[1]
      deList <- fitLinearModel(fit, conditions, es.design, deSize)
      message("Linear Models were fitted and saved in 'deList'.")
      
      a_deList <<- deList
      
      writeDifExprResultsToFiles(deList, conditions, dataSetSeries)
      message("Linear Models were written to files.")
      
      if (length(explanatoryTable) > 0)
        writeExplanatoryTableToFile(explanatoryTable, dataSetSeries)
      
      #deList <- readDifExprResultsFromFiles()
      #message("Linear Models were read from files and stored in 'deList'.")
      
      gseaResults <- geneSetEnrichmentAnalysis(deList, fileWithGenes)
      plots <<- gseaResults$gseaPlots
      gseaTableResults <<- gseaResults$gseaTableResults
      message("Gene set enrichment analysis was done.")
      
      if (length(gseaResults$gseaPlots) > 0) {
        writeGseaResults(gseaResults$gseaPlots, gseaResults$gseaTableResults,
                         conditions, dataSetSeries)
        message("Gene set enrichment analysis results were written to files.")
      }
    } else 
      message("Table isn't contained column with genes ID.")
    
  } else 
    message("Characteristics columns in the samples table don't exist or were unhelpful.
            So it need's to parse title column.")
} 
