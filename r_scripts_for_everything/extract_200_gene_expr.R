library(data.table)
library(dplyr)
source("./r_scripts_for_everything/extract_200_gene_handler.R")


gseTable <- read.table('some_documents/series_vs_gpl_tsv.csv', 
            header = T, sep = "\t", colClasses = c("character", "character"))

gmt <- list()

mineGSE <- readFile('some_documents/results_94.txt')
# mineGSE <- deleteRedundantFiles(mineGSE)
lastGSE <- ''

for (i in 1:length(mineGSE)) {
  message(i, ", ", mineGSE[i])
  
  tmp <- sub("results_94/", "", mineGSE[i])
  gse <- sub("/.*", '', tmp)

  if (!(grepl('-', gse)))
    gse <- getGSEwithGPL(gseTable, gse)
  
  pathGSE <- paste('dif_exprs/mm/', sub("/.*", '', tmp), sep='')
  pathFolder <- paste(pathGSE, '/dif_expression', sep='')
  # pathFolder <- paste('results_94/', sub("/.*", '', tmp), '/dif_expression_200', sep='')
  # pathFolder <- paste('results_94/', sub("/.*", '', tmp), '/dif_expression_400', sep='')
  # pathFolder <- paste('results_94/', sub("/.*", '', tmp), '/dif_expression_t_pvalue', sep='')
  # pathFolder <- paste('results_94/', sub("/.*", '', tmp), '/dif_expression_t_pvalue_log', sep='')
  
  if (length(gmt) == 0 || length(gmt) > 0 && !(isContainedGSE(gmt, gse))) {
    system(paste('mkdir ', pathGSE, sep=''), intern = T)
    system(paste('mkdir ', pathFolder, sep=''), intern = T)
  }
  
  if (grepl("explanatoryTable.tsv", mineGSE[i])) {
    system(paste('cp ', mineGSE[i], ' ', pathFolder, sep=''), intern = T)
  } else {
    condition <- sub(".*/", '', tmp)
    condition <- sub('.tsv', '', condition)
    condition <- strsplit(condition, split = "..vs..")[[1]]
    upCondition <- paste('up in \'', condition[1], '\' compared to \'', condition[2], '\'', sep = '')
    downCondition <- paste('up in \'', condition[2], '\' compared to \'', condition[1], '\'', sep = '')
    
    nameTables <- mineGSE[i]
    resultsGenes <- write200Genes(nameTables, pathFolder, upCondition, downCondition)
    # resultsGenes <- writeGenesDependsOnPvalueAndT(nameTables, pathFolder, upCondition, downCondition)
    
    if (!is.na(resultsGenes)) {
      if (!is.na(resultsGenes$universe) && (!is.na(resultsGenes$up) || !is.na(resultsGenes$down)) )
        if (i == 1) {
          gmt[[length(gmt)+1]] <- list("gse"=gse, "module"="universe", "genesID"=resultsGenes$universe)
          lastGSE <- gse
        } else if (gse != lastGSE) {
          gmt[[length(gmt)+1]] <- list("gse"=gse, "module"="universe", "genesID"=resultsGenes$universe)
          lastGSE <- gse
        }
        
      if (!is.na(resultsGenes$up))
        gmt[[length(gmt)+1]] <- list("gse"=gse, "module"=upCondition, "genesID"=resultsGenes$up)
      
      if (!is.na(resultsGenes$down))
        gmt[[length(gmt)+1]] <- list("gse"=gse, "module"=downCondition, "genesID"=resultsGenes$down)
    }
  }
}

writeGMT(gmt)
