downloadLibrariesAndSources <- function() {
  library(GEOquery)
  library(limma)
  library(gtools)
  library(stringr)
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


differentialExpression <- function(dataSetSeries, gpl) {

  gse <- getGEO(filename = dataSetSeries, getGPL = F)
  a_gse <<- gse
  dataSetSeries <- sub("(.*/)*", "", dataSetSeries) %>% 
    sub("_series_matrix.txt.gz", "", .)
  
  characteristics <- getCharacteristicsColumns(gse)
  conditionLists <- list()
  explanatoryTable <- data.table()
  
  if (length(characteristics) > 0) {
    message("There are characteristics columns in the samples table.")
    conStructure <- getConditionsFromCharacteristics(gse, characteristics)
    conditionLists <- conStructure$conditionsList
    explanatoryTable <- conStructure$explanatoryTable 
  }
  
  if (length(conditionLists) > 0 && length(unique(unlist(conditionLists))) > 1) {
    
    pData(gse)$condition <- fillGseConditionColumn(conditionLists)
    message("Column 'condition' was created and filled in the samples table.")
    
    gpl <- createGenesSymbolsTable(gpl)
    es <- collapseData(gse, gpl)
    
    a_gpl <<- gpl
    a_es <<- es
    
    message("Garbage was deleted from gene table.")

    if (!(is.na(max(exprs(es)))) && !(is.na(min(exprs(es))))) {
      
      if (max(exprs(es)) - min(exprs(es)) > 100)
        exprs(es) <- normalizeBetweenArrays(log2(exprs(es)+1), method="quantile")
      
      es.design <- model.matrix(~0+condition, data=pData(es))

      fit <- lmFit(es, es.design)
      conditions <- getConditionsForBuildingLinearModel(pData(gse)$condition)
      a_conditions <<- conditions
      a_des <<- es.design
      
      if (length(conditions) > 0) {
        message("Conditions combinations were received for filling of contrast matrix.")
    
        # Show received pairs of comparisons
        for (i in 1:length(conditions))
          message(conditions[[i]][1], " ", conditions[[i]][2])
        
        deSize <- dim(exprs(es))[1]
        
        deList <- fitLinearModel(fit, conditions, es.design, deSize)

        if (length(deList) > 0) {
          message("Linear Models were fitted and saved in 'deList'.")
          
          writeDifExprResultsToFiles(deList, conditions, dataSetSeries)
          message("Linear Models were written to files.")
          
          if (length(explanatoryTable) > 0)
            writeExplanatoryTableToFile(explanatoryTable, dataSetSeries)
    
        } else message("Linear Models can't be fitted.")
        
      } else message("There aren't conditions combinations with only one different condition.")
      
    } else message("Exprs table contains NA so linear Models can't be fitted.")
    
  } else 
    message("Characteristics columns in the samples table don't exist or were unhelpful.")
} 

