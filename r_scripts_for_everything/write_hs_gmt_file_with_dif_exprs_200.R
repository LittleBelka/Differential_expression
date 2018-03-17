library(data.table)
library(dplyr)
source("extract_200_gene_handler.R")


gseTable <- read.table('../some_documents/series_vs_gpl_tsv.csv', 
                       header = T, sep = "\t", colClasses = c("character", "character"))
gmt <- list()

mineGSE <- readFile('../some_documents/hs_dif_exprs_files_list.txt')

colTypes <- c("character", "character",
              "double", "double", "double", "double", "double", "double")

# for (i in 1:length(mineGSE)) {
for (i in 1:20) {
  message(i, ", ", mineGSE[i])
  
  newMineGSE <- paste('../', mineGSE[i], sep='')
  
  tmp <- sub("dif_exprs/hs/", "", mineGSE[i])
  gse <- sub("/.*", '', tmp)
  
  if (!(grepl('-', gse)))
    gse <- getGSEwithGPL(gseTable, gse)
  
  if (!grepl("explanatoryTable.tsv", mineGSE[i])) {
    tryCatch({
      t <- read.table(newMineGSE, header=T, sep="\t", 
                        fill = TRUE, quote = "", colClasses=colTypes)
      condition <- sub('.tsv', '', mineGSE[i]) %>% 
        sub('.*/', '', .)
      
      gmt[[length(gmt)+1]] <- list("gse"=gse, "module"=condition, "genesID"=t$rn)
      
    }, error = function(error_message) {
      message(error_message)
    })
  }
}

writeGMT(gmt)
