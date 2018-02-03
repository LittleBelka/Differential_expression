parseGPL <- function(gpl) {
  numberColumns <- list()
  for (i in 1:length(colnames(Table(gpl)))) {
    column = sapply(colnames(Table(gpl)[i]), tolower)[[1]]
    columnName = ''
    
    if (column == 'id')
      columnName = 'id'
    else if (column == 'gene symbol')
      columnName = 'gene_symbol'
    else if (column == 'gene id')
      columnName = 'gene_id'
    
    if (columnName != '')
      numberColumns[[length(numberColumns)+1]] <- list('number'=i, 'columnName'=columnName)
  }
  
  return(numberColumns)
}  


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


getListGPL <- function(listAllGPL, listHandledGPL) {
  listGPL <- c()
  tmpListHandledGPL <- c()
  
  for (i in 1:length(listHandledGPL)) {
    tmp <- sub('(.*/)*', '', listHandledGPL[i])
    tmp <- sub('.tsv', '', tmp)
    tmpListHandledGPL <- c(tmpListHandledGPL, tmp)
  }
  
  for (i in 1:length(listAllGPL)) {
    tmp <- sub('(.*/)*', '', listAllGPL[i])
    tmp <- sub('.annot.gz', '', tmp)
    if (!(tmp %in% tmpListHandledGPL)) listGPL <- c(listGPL, listAllGPL[i])
  }
  return(listGPL)
}









