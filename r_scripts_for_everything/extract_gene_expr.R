library(data.table)
library(dplyr)

source("./r_scripts_for_everything/extract_gene_handler.R")


gseTable <- read.table('some_documents/series_vs_gpl_tsv.csv', 
                       header = T, sep = "\t", colClasses = c("character", "character"))

gmt <- list()
gmtAnnotation <- data.table()

# mineGSE <- readFile('some_documents/results_94.txt')
mineGSE <- readFile('some_documents/rt.txt')
listGSEwithExplTables <- fillListMatchingGSEandExplTable(mineGSE)

lastGSE <- ''

moduleNumber <- 1

for (i in 1:length(mineGSE)) {
  message(i, ", ", mineGSE[i])
  
  # tmp <- sub("results_94/", "", mineGSE[i])
  tmp <- sub("rt/", "", mineGSE[i])
  gse <- sub("/.*", '', tmp)
  gseForExplTable <- gse
  gpl <- sub(".*-", '', gse)
  
  if (!(grepl('-', gse))) {
    resultGSE <- getGSEwithGPL(gseTable, gse)
    gse <- resultGSE$gseWithGPL
    gpl <- resultGSE$GPL
  } 
  
  # pathGSE <- paste('dif_exprs/mm/', sub("/.*", '', tmp), sep='')
  pathGSE <- paste('dif_exprs/rt/', sub("/.*", '', tmp), sep='')
  pathFolder <- paste(pathGSE, '/dif_expression', sep='')
  
  if (length(gmt) == 0 || length(gmt) > 0 && !(isContainedGSE(gmt, gse))) {
    system(paste('mkdir ', pathGSE, sep=''), intern = T)
    system(paste('mkdir ', pathFolder, sep=''), intern = T)
  }
  
  if (!grepl("explanatoryTable.tsv", mineGSE[i])) {
    
    conditions <- getUpDownConditions(gseForExplTable, tmp, listGSEwithExplTables)
    
    if (gse != lastGSE) moduleNumber <- 1
    
    
    resultsGenes <- write200Genes(mineGSE[i], pathFolder, moduleNumber)
    
    if (!is.na(resultsGenes)) {
      
      if (!is.na(resultsGenes$universe) && (!is.na(resultsGenes$up) || !is.na(resultsGenes$down)) ) {
        if (i == 1 || gse != lastGSE) {
          gmt[[length(gmt)+1]] <- list("gse"=gse, "module"=0, "genesID"=resultsGenes$universe)
          gmtAnnotation <- bind_rows(gmtAnnotation, 
                                     c(GSE=gse, GPL=gpl, Module_number=0, Module_name="universe"))
          lastGSE <- gse
        }
      }
      
      if (!is.na(resultsGenes$up)) {
        gmt[[length(gmt)+1]] <- list("gse"=gse, "module"=resultsGenes$moduleNumberUp, 
                                     "genesID"=resultsGenes$up)
        gmtAnnotation <- bind_rows(gmtAnnotation, c(GSE=gse, GPL=gpl, 
                                  Module_number=resultsGenes$moduleNumberUp, 
                                  Module_name=conditions$upCondition))
        moduleNumber <- moduleNumber + 1
      }
      
      if (!is.na(resultsGenes$down)) {
        gmt[[length(gmt)+1]] <- list("gse"=gse, "module"=resultsGenes$moduleNumberDown, 
                                     "genesID"=resultsGenes$down)
        gmtAnnotation <- bind_rows(gmtAnnotation, c(GSE=gse, GPL=gpl, 
                                  Module_number=resultsGenes$moduleNumberDown, 
                                  Module_name=conditions$downCondition))
        moduleNumber <- moduleNumber + 1
      }
      
      if (i == 4) {
        writeGMT(gmt, paste("rt.modules_with_numbers_", i, ".gmt", sep=''))
        writeGMTannotation(gmtAnnotation, paste("rt.modules_annotations_with_numbers_", i, ".tsv", sep=''))
      }
    }
  }
}

writeGMT(gmt, "rt.modules_with_numbers_last.gmt")
writeGMTannotation(gmtAnnotation, "rt.modules_annotations_with_numbers_last.tsv")
