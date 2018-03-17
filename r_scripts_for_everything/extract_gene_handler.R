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
  return(list("gseWithGPL"=gseWithGPL, "GPL"=result[[1]]))
}


fillListMatchingGSEandExplTable <- function(mineGSE) {
  listGSE <- list()
  for (i in 1:length(mineGSE)) {
    if (grepl("explanatoryTable.tsv", mineGSE[i])) {
      tmp <- sub("rt/", "", mineGSE[i])
      # tmp <- sub("results_94/", "", mineGSE[i])
      gse <- sub("/.*", '', tmp)
      
      table <- read.csv(mineGSE[i], sep = '\t', stringsAsFactors = F)
      listGSE[gse] <- list(table %>% 
        melt(id.vars = "marking") %>% 
        dcast(variable ~ marking) %>% 
        select(-variable))
    }
  }
  return(listGSE)
}


getUpDownConditions <- function(gse, gsePartOfPath, listGSEwithExplTables) {
  condition <- sub(".*/", '', gsePartOfPath)
  condition <- sub('.tsv', '', condition)
  condition <- strsplit(condition, split = "..vs..")[[1]]
  properConditions <- c()
  for (j in 1:length(condition)) {
    newCond <- paste('\'', condition[j], '\'', sep='') 
    if (gse %in% names(listGSEwithExplTables)) {
      newCond <- expandWithExplTable(paste('\'', condition[j], '\'', sep=''), 
                                     listGSEwithExplTables[gse][[1]])
    }
    properConditions <- c(properConditions, newCond)
  }
  
  upCondition <- paste('up in ', properConditions[1], ' compared to ', 
                       properConditions[2], sep = '')
  downCondition <- paste('up in ', properConditions[2], ' compared to ', 
                         properConditions[1], sep = '')
  
  return(list("upCondition"=upCondition, "downCondition"=downCondition))
}


expandWithExplTable <- function(condition, explTable) {
  
  for (m in names(explTable)) {
    parts <- strsplit(condition, m)[[1]]
    
    if (length(parts) > 1) {
      newCondition <- parts[1]
      
      for (i in 2:length(parts)) {
        partPrev <- strsplit(parts[i-1], '*')[[1]]
        partNext <- strsplit(parts[i], '*')[[1]]
        
        # message("PARTS: ", partPrev, " !! ", partNext, " !m: ", m)
        
        if (length(partPrev) != 0 && length(partNext) != 0 &&
            (partPrev[length(partPrev)] == "'" && partNext[1] == "'" ||
            partPrev[length(partPrev)] == "_" && partNext[1] == "'" ||
            partPrev[length(partPrev)] == "'" && partNext[1] == "_" ||
            partPrev[length(partPrev)] == "_" && partNext[1] == "_")) {
          
          newCondition <- paste(newCondition, gsub(' ', '.', explTable[m]), sep='')
        } else {
          newCondition <- paste(newCondition, m, sep='')
        }
        newCondition <- paste(newCondition, parts[i], sep='')
      }
      
      newCondition_a <<- newCondition
      condition <- newCondition
      # message("NEW_CONDITION: ", newCondition, " !condition: ", condition)
    }
  }
  return(condition) 
}


write200Genes <- function(nameTables, pathFolder, moduleNumber) {
  colTypes <- c("character", "character",
                "double", "double", "double", "double", "double", "double")
  tryCatch({
    gse <- read.table(nameTables, header=T, sep="\t", 
                      fill = TRUE, quote = "", colClasses=colTypes)
    gse <- subset(gse, !is.na(gse$symbol)) 
    
    up <- gse[order(gse$t, decreasing=TRUE),][1:200,]
    down <- gse[order(gse$t),][1:200,]
    
    moduleNumberUp <- -1
    moduleNumberDown <- -1
    
    result <- writeTSV(up, pathFolder, moduleNumber)
    
    if (result == "success") {
      moduleNumberUp <- moduleNumber
      moduleNumber <- moduleNumber + 1
    }
    
    result <- writeTSV(down, pathFolder, moduleNumber)
    
    if (result == "success") moduleNumberDown <- moduleNumber
    
    return(list("up"=up$rn, "down"=down$rn, "universe"=gse$rn, 
                "moduleNumberUp"=moduleNumberUp, "moduleNumberDown"=moduleNumberDown))
  }, error = function(error_message) {
    message(error_message)
  })
  return(NA)
}


writeTSV <- function(gse, pathFolder, moduleNumber) {
  if (nrow(gse) > 0) {
    fileName <- paste(pathFolder, '/', moduleNumber, '.tsv', sep = '')
    fwrite(gse, fileName, sep = '\t')
    return("success")
  }
  return("fail")
}


writeGMT <- function(gmt, fileName) {
  if (length(gmt) > 0) {
    fileConn<-file(paste("some_documents/", fileName, sep=''))
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

writeGMTannotation <- function(gmtAnnotation, fileName) {
  write.table(gmtAnnotation, paste("some_documents/", fileName, sep=''), 
              sep="\t", row.names = F)
}


