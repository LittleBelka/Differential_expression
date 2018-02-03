library(GEOquery)
library(data.table)
# for handle of reading gpl
source("./r_scripts_for_everything/read_gpl_handler.R")

listAllGPL <- readFile("r_scripts_for_everything/list_all_gpl.txt")
listHandledGPL <- readFile("r_scripts_for_everything/list_current_gpl.txt")

listGPL <- getListGPL(listAllGPL, listHandledGPL)

listGPL <- system('find data -type f', intern = T)
# listGPL <- system('find ../geo/platforms/ -type f', intern = T)
#listGPL <- list("data/GPL10666.annot.gz")

for (k in 1:length(listGPL)) {

  gpl <- getGEO(filename = listGPL[[k]])
  numberColumns <- parseGPL(gpl)
  numbers <- c(sapply(numberColumns, function(x) x$number))
  names <- c(sapply(numberColumns, function(x) x$columnName))
  
  gplTable <- data.frame(Table(gpl)[numbers], stringsAsFactors = F)
  colnames(gplTable) <- names
  
  gplTable <- subset(gplTable, gplTable$gene_symbol != '' & gplTable$gene_id != '')
  tmpTable <- subset(gplTable, grepl("///", gplTable$gene_symbol))
  gplTable <- subset(gplTable, !grepl("///", gplTable$gene_symbol))
  gplTable <- subset(gplTable, !is.na(gplTable$gene_symbol))
  
  gplTable$id <- as.character(gplTable$id)
  gplTable$gene_symbol <- as.character(gplTable$gene_symbol)
  gplTable$gene_id <- as.character(gplTable$gene_id)
  
  if (nrow(tmpTable) > 0) {
    tmpTable$id <- as.character(tmpTable$id)
    tmpTable$gene_symbol <- as.character(tmpTable$gene_symbol)
    tmpTable$gene_id <- as.character(tmpTable$gene_id)
    
    for (i in 1:nrow(tmpTable)) {
      geneSymbol <- strsplit(as.character(tmpTable[i, 2]), '///')
      geneID <- strsplit(as.character(tmpTable[i, 3]), '///')
      
      if (length(geneSymbol[[1]]) > 0 && length(geneID[[1]]) > 0 && 
          length(geneSymbol[[1]]) == length(geneID[[1]])) {
        
        for (j in 1:length(geneSymbol[[1]])) {
          # e <- c(tmpTable[i, 1], geneSymbol[[1]][j], geneID[[1]][j])
          # gplTable <- rbind(gplTable, e, stringsAsFactors=F)
          #gplTable[nrow(gplTable) + 1,] <- c(tmpTable[i, 1], geneSymbol[[1]][j], geneID[[1]][j])
          #gplTable[nrow(gplTable) + 1,] <- c(NA, NA, NA)
          gplTable[nrow(gplTable) + 1,] <- c(as.character(tmpTable[i, 1]), 
                                             as.character(geneSymbol[[1]][j]), as.character(geneID[[1]][j]))
        }
      }
    }
  }
  
  gplName <- sub('(.*/)*', '', listGPL[[k]])
  gplName <- sub('.annot.gz', '', gplName)
  
  fwrite(gplTable, paste('../gpl/', gplName, '.tsv', sep = ''), sep = '\t')
}
warnings()
print(tail(gplTable))
