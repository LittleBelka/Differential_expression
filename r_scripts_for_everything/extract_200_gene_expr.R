library(data.table)
library(dplyr)
source("./r_scripts_for_everything/extract_200_gene_handler.R")


gseTable <- read.table('some_documents/series_vs_gpl_tsv.csv', 
            header = T, sep = "\t", colClasses = c("character", "character"))

gmt <- list()

mineGSE <- readFile('some_documents/results_93.txt')
mineGSE <- deleteRedundantFiles(mineGSE)

for (i in 1:length(mineGSE)) {
  message(i, ", ", mineGSE[i])
  
  tmp <- sub("results_158/", "", mineGSE[i])
  gse <- sub("/.*", '', tmp)

  if (!(grepl('-', gse)))
    gse <- getGSEwithGPL(gseTable, gse)
  
  pathFolder <- paste('results_158/', sub("/.*", '', tmp), '/dif_expression_200', sep='')
  
  if (length(gmt) == 0 || length(gmt) > 0 && !(isContainedGSE(gmt, gse))) 
    system(paste('mkdir ', pathFolder, sep=''), intern = T)
  
  condition <- sub(".*/", '', tmp)
  condition <- sub('.tsv', '', condition)
  
  nameTables <- mineGSE[i]
  resultsGenes <- write200Genes(nameTables, pathFolder, condition)
  
  if (!is.na(resultsGenes)) {
    gmt[[length(gmt)+1]] <- list("gse"=gse, "module"=paste(condition, '.up', sep = ''), 
                                 "genesID"=resultsGenes$up)
    
    gmt[[length(gmt)+1]] <- list("gse"=gse, "module"=paste(condition, '.down', sep = ''), 
                                 "genesID"=resultsGenes$down)
  }
}

writeGMT(gmt)
