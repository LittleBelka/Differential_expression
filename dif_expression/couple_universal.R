downloadLibrariesAndSources <- function() {
  library(GEOquery)
  library(limma)
  # library(org.Mm.eg.db)
  # library(org.Rn.eg.db)
  # library(org.Hs.eg.db)
  library(gtools)
  library(stringr)
  # library(fgsea)
  library(ggplot2)
  library(data.table)
  library(dplyr)
  library(magrittr)
  library(tibble)
  # for makeContrasts2
  source("https://raw.githubusercontent.com/assaron/r-utils/master/R/exprs.R")
  # for conditions handling
  source("./dif_expression/couple_universal_handler.R")
  # for gene set enrichment analysis
  source("./dif_expression/gsea.R")
}


differentialExpression <- function(dataSetSeries, gpl, fileWithGenes) {

  gse <- getGEO(filename = dataSetSeries, getGPL = F)
  
  dataSetSeries <- sub("(.*/)*", "", dataSetSeries) %>% 
    sub("_series_matrix.txt.gz", "", .)
  
  characteristics <- getCharacteristicsColumns(gse)
  conditionLists <- list()
  explanatoryTable <- data.table()
  
  if (length(characteristics) > 0) {
    message("There are characteristics columns in the samples table.")
    conStructure <- getConditionsFromCharacteristics(gse, characteristics)
    a_conStructure <<- conStructure
    conditionLists <- conStructure$conditionsList
    explanatoryTable <- conStructure$explanatoryTable 
  }
  
  if (length(conditionLists) > 0 && length(unique(unlist(conditionLists))) > 1) {
    pData(gse)$condition <- fillGseConditionColumn(conditionLists)
    message("Column 'condition' was created and filled in the samples table.")
    
    a_gse <<- gse
    a_gpl <<- gpl
    gpl <- createGenesSymbolsTable(gpl)
    es <- collapseData(gse, gpl)
    
    a_es <<- es
    message("Garbage was deleted from gene table.")
    # fData(es) <- data.frame(row.names = rownames(es))
    # database <- getDatabaseForMapping(gse)
    # fData(es)$symbol <- mapIds(database, keys = rownames(es), column = "SYMBOL", keytype = "ENTREZID")

    if (max(exprs(es)) - min(exprs(es)) > 100)
      exprs(es) <- normalizeBetweenArrays(log2(exprs(es)+1), method="quantile")
    
    es.design <- model.matrix(~0+condition, data=pData(es))

    fit <- lmFit(es, es.design)
    conditions <- getConditionsForBuildingLinearModel(pData(gse)$condition)

    message("Conditions combinations were received for filling of contrast matrix.")

    # Show received pairs of comparisons
    for (i in 1:length(conditions))
      message(conditions[[i]][1], " ", conditions[[i]][2])
    
    deSize <- dim(exprs(es))[1]
    
    a_conditions <<- conditions
    a_fit <<- fit
    a_es.design <<- es.design
    deSize <<- deSize
    
    deList <- fitLinearModel(fit, conditions, es.design, deSize)
    if (length(deList) > 0) {
      message("Linear Models were fitted and saved in 'deList'.")
      
      writeDifExprResultsToFiles(deList, conditions, dataSetSeries)
      message("Linear Models were written to files.")
      
      if (length(explanatoryTable) > 0)
        writeExplanatoryTableToFile(explanatoryTable, dataSetSeries)
  
      #deList <- readDifExprResultsFromFiles()
      #message("Linear Models were read from files and stored in 'deList'.")
  
      # gseaResults <- geneSetEnrichmentAnalysis(deList, fileWithGenes)
      # plots <<- gseaResults$gseaPlots
      # gseaTableResults <<- gseaResults$gseaTableResults
      # message("Gene set enrichment analysis was done.")
      # 
      # if (length(gseaResults$gseaPlots) > 0) {
      #   writeGseaResults(gseaResults$gseaPlots, gseaResults$gseaTableResults,
      #                    conditions, dataSetSeries)
      #   message("Gene set enrichment analysis results were written to files.")
      # }
    } else message("Linear Models can't be fitted.")
    
  } else 
    message("Characteristics columns in the samples table don't exist or were unhelpful.
            So it need's to parse title column.")
} 


findGSEvsGPL <- function() {
  table <- read.csv('some_documents/series_vs_gpl_tsv.csv', 
                                      sep = '\t', stringsAsFactors = F)
  result <- table %>% 
    filter(GPL != '') %>% 
    sample_n(3) %>% 
    select(Series_path, GPL) %>% 
    rename(dataSetSeries = Series_path, availableGPL = GPL)
  
  return(result)
}

