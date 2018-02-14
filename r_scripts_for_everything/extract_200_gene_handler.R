readFile <- function(fileName) {
  con = file(fileName, "r")
  strings <- c()
  while (TRUE) {
    line = readLines(con, n = 1)
    strings <- c(strings, line) 
    if (length(line) == 0) break
  }
  close(con)
  return(strings)
}


isContainedGSE <- function(gmt, gse) {
  for(i in 1:length(gmt)) { 
    if (gse %in% gmt[[i]]$gse) return(T)
  }
  return(F)
}


getGSEwithGPL <- function(gseTable, gse) {
  result <- gseTable %>% 
    filter(GPL != '') %>% 
    filter(Series == gse) %>% 
    select(Platforms)
  gseWithGPL <- paste(gse, result[[1]], sep='-')
}


deleteRedundantFiles <- function(mineGSE) {
  gse <- c()
  for (g in mineGSE) {
    if (!(grepl("explanatoryTable.tsv", g))) gse <- c(gse, g)
  }
  return(gse)
}


write200Genes <- function(nameTables, pathFolder, condition) {
  colTypes <- c("character", "character",
                "double", "double", "double", "double", "double", "double")
  tryCatch({
    gse <- read.table(nameTables, header=T, sep="\t", colClasses=colTypes)
    gse <- subset(gse, !is.na(gse$symbol))
    up <- gse[order(gse$t, decreasing=TRUE),][1:200,]
    down <- gse[order(gse$t),][1:200,]
    
    writeCSV(up, pathFolder, condition, '.up')
    writeCSV(down, pathFolder, condition, '.down')
    
    return(list("up"=up$rn, "down"=down$rn))
  }, error = function(error_message) {
    message(error_message)
  })
  return(NA)
}


writeCSV <- function(gse, pathFolder, condition, regulation) {
  
  fileName <- paste(pathFolder, '/', condition, regulation, '.tsv', sep = '')
  fwrite(gse, fileName, sep = '\t')
}


writeGMT <- function(gmt) {
  fileConn<-file("some_documents/modules.gmt")
  toWrite <- c()
  
  for (i in 1:length(gmt)) {
    id <- paste(gmt[[i]]$genesID, collapse=",")
    s <- paste(gmt[[i]]$gse, "#", gmt[[i]]$module, 
               "    ", id, sep='')
    toWrite <- c(toWrite, s)
  }
  
  writeLines(toWrite, fileConn)
  close(fileConn)
}




