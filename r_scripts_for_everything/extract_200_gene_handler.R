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
  return(gseWithGPL)
}


deleteRedundantFiles <- function(mineGSE) {
  gse <- c()
  for (g in mineGSE) {
    if (!(grepl("explanatoryTable.tsv", g))) gse <- c(gse, g)
  }
  return(gse)
}


write200Genes <- function(nameTables, pathFolder, upCondition, downCondition) {
  colTypes <- c("character", "character",
                "double", "double", "double", "double", "double", "double")
  tryCatch({
    gse <- read.table(nameTables, header=T, sep="\t", 
                      fill = TRUE, quote = "", colClasses=colTypes)
    gse <- subset(gse, !is.na(gse$symbol)) 
    gse <- gse[order(gse$AveExpr, decreasing=TRUE),][1:7000,]

    # up <- gse[order(gse$t, decreasing=TRUE),][1:200,]
    # down <- gse[order(gse$t),][1:200,]
    up <- gse[order(gse$t, decreasing=TRUE),][1:400,]
    down <- gse[order(gse$t),][1:400,]
    
    writeCSV(up, pathFolder, upCondition)
    writeCSV(down, pathFolder, downCondition)
    
    return(list("up"=up$rn, "down"=down$rn))
  }, error = function(error_message) {
    message(error_message)
  })
  return(NA)
}


writeGenesDependsOnPvalueAndT <- function(nameTables, pathFolder, upCondition, downCondition) {
  colTypes <- c("character", "character",
                "double", "double", "double", "double", "double", "double")
  tryCatch({
    gse <- read.table(nameTables, header=T, sep="\t", 
                      fill = TRUE, quote = "", colClasses=colTypes)
    gse <- subset(gse, !is.na(gse$symbol)) 
    gse <- subset(gse, !is.na(gse$t))
    gse <- subset(gse, !is.na(gse$adj.P.Val))
    gse <- gse[order(gse$AveExpr, decreasing=TRUE),][1:7000,]
    
    gse <- subset(gse, abs(gse$logFC) > 0.5)
    
    up <- subset(gse, gse$t > 0 & gse$adj.P.Val < 0.05)
    down <- subset(gse, gse$t < 0 & gse$adj.P.Val < 0.05)
    
    writeCSV(up, pathFolder, upCondition)
    writeCSV(down, pathFolder, downCondition)
    
    if (nrow(up) > 0 && nrow(down) > 0) {
      return(list("up"=up$rn, "down"=down$rn))
    } else if (nrow(up) > 0 && nrow(down) == 0) {
      return(list("up"=up$rn, "down"=NA))
    } else if (nrow(up) == 0 && nrow(down) > 0) {
      return(list("up"=NA, "down"=down$rn))
    }
    
  }, error = function(error_message) {
    message(error_message)
  })
  return(NA)
}


writeCSV <- function(gse, pathFolder, condition) {
  if (nrow(gse) > 0) {
    fileName <- paste(pathFolder, '/', condition, '.tsv', sep = '')
    fwrite(gse, fileName, sep = '\t')
  }
}


writeGMT <- function(gmt) {
  if (length(gmt) > 0) {
    fileConn<-file("some_documents/modules_t_pvalue_log.gmt")
    # fileConn<-file("some_documents/modules_400_genes.gmt")
    toWrite <- c()
    
    for (i in 1:length(gmt)) {
      id <- paste(gmt[[i]]$genesID, collapse=",")
      s <- paste(gmt[[i]]$gse, '#', gmt[[i]]$module, 
                 '\t' , id, sep='')
      toWrite <- c(toWrite, s)
    }
    
    writeLines(toWrite, fileConn)
    close(fileConn)
  }
}




