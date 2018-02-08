getCharacteristicsColumns <- function(gse) {
  col <- colnames(pData(gse))
  characteristics <- c()
  for (ch in col) {
    if (grepl("characteristics", ch)) characteristics <- c(characteristics, ch)
  }
  return(characteristics)
}


getConditionsFromTitle <- function() {}


createMarkersForComplicatedConditions <- function() {
  usedMarkers <- list()
  for (i in 1:26) usedMarkers[intToUtf8(64+i)] <- 0
  
  extraMarkers <- permutations(
    n=10, r=2, v=names(usedMarkers)[1:10], set=T, repeats.allowed=T)
  
  for (i in 1:(length(extraMarkers))/2) {
    cond <- paste(extraMarkers[[i,1]], extraMarkers[[i,2]], sep="")
    usedMarkers[cond] <- 0
  }
  return(usedMarkers)
}


getConditionsFromCharacteristics <- function(gse, characteristics) {
  uniqList <- list()
  conditionsList <- list()
  explanatoryTable <- data.table()
  usedMarkers <- createMarkersForComplicatedConditions()
  
  for (ch in characteristics) {
    uniqChar <- unique(pData(gse)[ch])
    
    if (length(rownames(uniqChar)) > 1) {
      conStructure <- getConditionsFromUniqValues(pData(gse)[ch], conditionsList, 
                                                  explanatoryTable, usedMarkers)
      conditionsList <- conStructure$conditionsList
      explanatoryTable <- conStructure$explanatoryTable
      usedMarkers <- conStructure$usedMarkers
    }
  }
  return(list("conditionsList"=conditionsList, "explanatoryTable"=explanatoryTable))
}


getConditionsFromUniqValues <- function(charColumns, conditionsList, 
                                              explanatoryTable, usedMarkers) {
  values <- sub("^.*: ", "", charColumns[[1]])
  numerics <- is.na(suppressWarnings(as.numeric(values)))
  
  if (length(numerics[numerics==FALSE]) == 0) {
    tmpCondition <- c()
    tmpComplicatedCondition <- list()
    
    for (v in values) {
      oldV <- v
      complicatedCondition <- FALSE
      
      if (grepl("\\(.*\\)", v)) v <- tryGetConditionFromBrackets(v)
      
      if (isTooComplicatedCondition(v)) complicatedCondition <- TRUE
      
      if (complicatedCondition) {
        resultCondition <- handleComplicatedCondition(tmpComplicatedCondition, 
                                          oldV, explanatoryTable, usedMarkers)
        v <- resultCondition$condition
        tmpComplicatedCondition <- resultCondition$tmpComplicatedCondition
        explanatoryTable <- resultCondition$explanatoryTable
        usedMarkers <- resultCondition$usedMarkers
      } else 
        v <- gsub("[^a-zA-Z0-9]*", '\\1', v) 
      
      if (grepl("[a-zA-Z]+", v)) tmpCondition <- c(tmpCondition, v)
    }
    conditionsList[[length(conditionsList)+1]] <- tmpCondition
  }
  return(list("conditionsList"=conditionsList, "explanatoryTable"=explanatoryTable,
              "usedMarkers"=usedMarkers))
}


tryGetConditionFromBrackets <- function(condition) {
  newCondition <- str_replace(condition, '.*\\((.*)\\).*', '\\1')
  tmp <- gsub("[^a-zA-Z]*", '\\1', newCondition)
  if (nchar(tmp) > 0)
    return(newCondition)
  return(condition)
}


handleComplicatedCondition <- function(tmpComplicatedCondition, cond, 
                                              explanatoryTable, usedMarkers) {
  cond <- gsub("[^a-zA-Z0-9 \\+-]*", '\\1', cond)
  newConditionName <- isContainedInConditionList(tmpComplicatedCondition, cond)
  if (newConditionName != "") {
    newCond <- newConditionName 
  } else {
    markerUse <- getFreeMarker(usedMarkers)
    tmpComplicatedCondition[cond] <- markerUse$marker
    usedMarkers <- markerUse$usedMarkers
    
    explanatoryTable <- rbindlist(list(explanatoryTable, data.table(
            "characteristics"=cond, "marking"=tmpComplicatedCondition[cond])))
  }
  return(list("condition"=tmpComplicatedCondition[cond], 
              "tmpComplicatedCondition"=tmpComplicatedCondition, 
              "explanatoryTable"=explanatoryTable,
              "usedMarkers"= usedMarkers))
}


isContainedInConditionList <- function(tmpComplicatedCondition, cond) {
  if (cond %in% names(tmpComplicatedCondition)) 
    return(tmpComplicatedCondition[cond])
  
  return("")
}


getFreeMarker <- function(usedMarkers) {
  i = 0
  found <- FALSE
  marker <- ""
  while(!found && length(usedMarkers) > i) {
    i = i + 1
    if (usedMarkers[i] == 0) {
      found <- TRUE
      marker <- names(usedMarkers)[i]
      usedMarkers[i] <- 1
    }
  }
  return(list("usedMarkers"=usedMarkers, "marker"=marker))
}


isTooComplicatedCondition <- function(condition) {
  if (length(unlist(gregexpr(pattern = "\\+", condition))) > 1
      || length(unlist(gregexpr(pattern = "-", condition))) > 1
      || length(unlist(gregexpr(pattern = " ", condition))) > 1)
    return(TRUE)
  return(FALSE)
}


fillGseConditionColumn <- function(conditionList) {
  conditions <- c()
  conditionList <- matrix(unlist(conditionList), 
                          ncol = length(conditionList[[1]]), byrow = TRUE)
  conditionList <- split(conditionList, c(col(conditionList)))
  
  for (i in 1:length(conditionList)) {
    tmpCondition <- ""
    
    for (j in 1:length(conditionList[[i]]))
      tmpCondition <- paste(tmpCondition, conditionList[[i]][j], sep="_")
    
    tmpCondition <- str_replace(tmpCondition, "_$", "")
    conditions <- c(conditions, str_replace(tmpCondition, "^_", ""))
  }
  return(conditions)
}


getConditionsForBuildingLinearModel <- function(conditions) {
  tmpConditions <- unique(conditions)
  conditions <- list()
  tmpConditions <- combinations(
          n=length(tmpConditions), r=2, v=tmpConditions, set=T, repeats.allowed=F)
  
  for (i in 1:(length(tmpConditions)/2)) {
    firstCon <- strsplit(tmpConditions[[i,1]], "_")
    secondCon <- strsplit(tmpConditions[[i,2]], "_")
    k <- 0
    for (j in 1:length(firstCon[[1]])) {
      if (firstCon[[1]][j] != secondCon[[1]][j]) k <- k + 1
    }
    if (k == 1) conditions[[length(conditions)+1]] <- 
      list("firstCon"=tmpConditions[[i,1]], "secondCon"=tmpConditions[[i,2]])
  }
  return(conditions)
}


# provideValidOfSomeColumns <- function(gse) {
#   successValidation <- FALSE
#   if ("ENTREZ_GENE_ID" %in% names(fData(gse))) {
#     successValidation <- TRUE
#   } else {
#     columnNames <- names(fData(gse))
#     
#     for (i in 1:length(columnNames)) {
#       if (sapply(columnNames[i], tolower)[[1]] == "entrez_gene_id") {
#         names(fData(gse))[names(fData(gse)) == columnNames[i]] <- "ENTREZ_GENE_ID"
#         successValidation <- TRUE
#       } else if (sapply(columnNames[i], tolower)[[1]] == "gene") {
#         names(fData(gse))[names(fData(gse)) == columnNames[i]] <- "GENE"
#         numerics <- is.na(suppressWarnings(as.numeric(fData(gse)$GENE)))
#         
#         if (sum(!is.na(fData(gse)$GENE)) == length(numerics[numerics==FALSE]))
#           successValidation <- TRUE
#       }
#     }
#   }
#   return(list("fData"=fData(gse), "successValidation"=successValidation))
# }


# getDatabaseForMapping <- function(gse) {
#   col <- colnames(pData(gse))
#   curOrganism <- ""
#   for (ch in col) 
#     if (grepl("organism", ch)) curOrganism <- tolower(pData(a_gse)[ch][[1]][1])
#   
#   database <- org.Hs.eg.db
#   
#   if (curOrganism != "") {
#     if (curOrganism == "mus musculus") database <- org.Mm.eg.db
#     if (curOrganism == "rattus norvegicus") database <- org.Rn.eg.db
#   }
#   
#   return(database)
# }


fitLinearModel <- function(fit, conditions, design, deSize) {
  deList <- list()
  for (i in 1:length(conditions)) {
    
    contrasts <- makeContrasts2(
               c("condition", conditions[[i]]$firstCon, conditions[[i]]$secondCon), 
               levels=design)
    
    fit2 <- contrasts.fit(fit, contrasts)
    
    df.residual <- unique(fit2$df.residual)
    if ((length(df.residual) == 1 && df.residual[1] != 0) ||
        (length(df.residual) != 1)) {
      
      fit2 <- eBayes(fit2)
      
      deList[[i]] <- data.table(
        topTable(fit2, adjust.method="BH", number=deSize, sort.by = "B"), 
        keep.rownames = T)
    }
  }
  return(deList)
}


writeExplanatoryTableToFile <- function(explanatoryTable, dataSetSeries) {
  nameFile <- paste("./results/", dataSetSeries, 
                    "/dif_expression/explanatoryTable.tsv", sep="")
  fwrite(explanatoryTable, nameFile, sep = "\t")
}


writeDifExprResultsToFiles <- function(deList, conditions, dataSetSeries) {
  dir.create(file.path("./results/", dataSetSeries), showWarnings = FALSE)
  filePath <- paste("./results/", dataSetSeries, "/", sep="")
  dir.create(file.path(filePath, "dif_expression"), showWarnings = FALSE)
  
  for (i in 1:length(conditions)) {
    nameFile <- paste(filePath,
                      "/dif_expression/",
                      conditions[[i]]$firstCon, ".vs.", 
                      conditions[[i]]$secondCon, ".tsv", sep="")
    write.table(deList[[i]], nameFile, sep="\t", quote = F, row.names = F)
  }
}


# readDifExprResultsFromFiles <- function() {
#   colTypes <- c("character", "character", 
#                 "double", "double", "double", "double", "double","double")
#   files <- dir(path = "./results/dif_expression", 
#                full.names = TRUE, recursive = TRUE)
#   deList <- list()
#   
#   for (i in 1:length(files))
#     deList[[i]] <- read.table(files[i], header=T, colClasses=colTypes)
#   
#   return(deList)
# }


createGenesSymbolsTable <- function(gpl) {
  gpl <- read.table(gpl, header = T, sep = "\t",
                    colClasses = c("character", "character", "integer"))
  gpl <- gpl[!duplicated(gpl$id),]
  gpl <- gpl[!duplicated(gpl$gene_id),]
  
  gpl <- gpl %>% 
    set_rownames(.$id) %>% 
    select(ID = id, ENTREZ_GENE_ID = gene_id, symbol = gene_symbol)

  return(gpl)
}


collapseData <- function(gse, gpl, FUN=median) {
  ranks <- apply(exprs(gse), 1, FUN)
  ranks <- data.frame(r = ranks, i = seq_along(ranks))
  table <- inner_join(rownames_to_column(gpl), rownames_to_column(ranks), 
           by="rowname") %>% 
           mutate(j = seq_along(symbol))
  t <- table[order(table$r, decreasing=T), ]
  keep <- t$i
  res <- gse[keep, ]
  rownames(res) <- table$ENTREZ_GENE_ID[t$j]
  return(res)
}
