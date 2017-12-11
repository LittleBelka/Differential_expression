getCharacteristicsColumns <- function(gse) {
  col <- colnames(pData(gse))
  characteristics <- c()
  for (ch in col) {
    if (grepl("characteristics", ch)) characteristics <- c(characteristics, ch)
  }
  return(characteristics)
}


getConditionsFromTitle <- function() {}


getConditionsFromCharacteristics <- function(gse, characteristics) {
  uniqList <- list()
  conditionsList <- list()
  explanatoryTable <- data.table()
  
  for (ch in characteristics) {
    uniqChar <- unique(pData(gse)[ch])
    
    if (length(rownames(uniqChar)) > 1) {
      conStructure <- getConditionsFromUniqValues(gse, pData(gse)[ch], 
                                      conditionsList, uniqChar, explanatoryTable)
      conditionsList <- conStructure$conditionsList
      explanatoryTable <- conStructure$explanatoryTable
    }
  }
  return(list("conditionsList"=conditionsList, "explanatoryTable"=explanatoryTable))
}


getConditionsFromUniqValues <- function(gse, charColumns, conditionsList, uniqChar, 
                                                                explanatoryTable) {
  values <- sub("^.*: ", "", charColumns[[1]])
  numerics <- is.na(suppressWarnings(as.numeric(values)))
  
  if (length(numerics[numerics==FALSE]) == 0) {
    tmpCondition <- c()
    tmpComplicatedCondition <- createTmpComplicatedConditionList(uniqChar)
    
    for (v in values) {
      oldV <- v
      complicatedCondition <- FALSE
      
      if (grepl("\\(.*\\)", v)) v <- tryGetConditionFromBrackets(v)
      
      if (isTooComplicatedCondition(v)) complicatedCondition <- TRUE
      
      if (complicatedCondition) {
        resultCondition <- handleComplicatedCondition(tmpComplicatedCondition, 
                                                      oldV, explanatoryTable)
        v <- resultCondition$condition
        tmpComplicatedCondition <- resultCondition$tmpComplicatedCondition
        explanatoryTable <- resultCondition$explanatoryTable
      } else 
        v <- gsub("[^a-zA-Z0-9]*", '\\1', v) 
      
      if (grepl("[a-zA-Z]+", v)) tmpCondition <- c(tmpCondition, v)
    }
    conditionsList[[length(conditionsList)+1]] <- tmpCondition
  }
  return(list("conditionsList"=conditionsList, "explanatoryTable"=explanatoryTable))
}


tryGetConditionFromBrackets <- function(condition) {
  newCondition <- str_replace(condition, '.*\\((.*)\\).*', '\\1')
  tmp <- gsub("[^a-zA-Z]*", '\\1', newCondition)
  if (nchar(tmp) > 0)
    return(newCondition)
  return(condition)
}


createTmpComplicatedConditionList <- function(uniqChar) {
  tmpComplicatedCondition <- list()
  lenUniqChar <- length(rownames(uniqChar))
  letters <- c()
  for (i in 1:26) letters <- c(letters, intToUtf8(64+i))
  i = 1
  
  while (lenUniqChar >= i && i <= 26) {
    tmpComplicatedCondition[[i]] <- list("condition"="", "replacement"=letters[i])
    i = i + 1
  }
  if (lenUniqChar > 26) 
    tmpComplicatedCondition <- condListWithPermutations(
                                tmpComplicatedCondition, letters, lenUniqChar, i)
  
  return(tmpComplicatedCondition)
}


condListWithPermutations <- function(tmpComplicatedCondition, letters, 
                                                            lenUniqChar, i) {
  tmpConditions <- permutations(
    n=length(letters[1:6]), r=2, v=letters[1:6], set=T, repeats.allowed=T)
  
  while (lenUniqChar >= i) {
    cond <- paste(tmpConditions[[i,1]], tmpConditions[[i,2]], sep="")
    tmpComplicatedCondition[[i]] <- list("condition"="", 
                                         "replacement"=cond)
    i = i + 1
  }
  return(tmpComplicatedCondition)
}


handleComplicatedCondition <- function(tmpComplicatedCondition, cond, 
                                                        explanatoryTable) {
  cond <- gsub("[^a-zA-Z0-9 \\+-]*", '\\1', cond) 
  newConditionName <- isContainedInConditionList(tmpComplicatedCondition, cond)
  if (newConditionName != "") {
    newCond <- newConditionName 
  } else {
    newConditionIndex <- getFreeListPosition(tmpComplicatedCondition)
    tmpComplicatedCondition[[newConditionIndex]]$condition <- cond
    newCond <- tmpComplicatedCondition[[newConditionIndex]]$replacement
    
    explanatoryTable <- rbindlist(list(explanatoryTable, data.table(
                                  "characteristics"=cond, "marking"=newCond)))
  }
  return(list("condition"=newCond, 
              "tmpComplicatedCondition"=tmpComplicatedCondition, 
              "explanatoryTable"=explanatoryTable))
}


getFreeListPosition <- function(tmpComplicatedCondition) {
  i = 0
  found <- FALSE
  while(!found && length(tmpComplicatedCondition) > i) {
    i = i + 1
    if (tmpComplicatedCondition[[i]]$condition == "") found <- TRUE
  }
  return(i)
}


isContainedInConditionList <- function(tmpComplicatedCondition, condition) {
  for (i in 1:length(tmpComplicatedCondition)) {
    if (length(tmpComplicatedCondition[[i]]$condition) > 0 
        && tmpComplicatedCondition[[i]]$condition == condition) 
      return(tmpComplicatedCondition[[i]]$replacement)
  }
  return("")
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


provideValidOfSomeColumns <- function(gse) {
  successValidation <- FALSE
  if ("ENTREZ_GENE_ID" %in% names(fData(gse))) {
    successValidation <- TRUE
  } else {
    columnNames <- names(fData(gse))
    
    for (i in 1:length(columnNames)) {
      if (sapply(columnNames[i], tolower)[[1]] == "entrez_gene_id") {
        names(fData(gse))[names(fData(gse)) == columnNames[i]] <- "ENTREZ_GENE_ID"
        successValidation <- TRUE
      } else if (sapply(columnNames[i], tolower)[[1]] == "gene") {
        names(fData(gse))[names(fData(gse)) == columnNames[i]] <- "GENE"
        numerics <- is.na(suppressWarnings(as.numeric(fData(a_gse)$GENE)))
        
        if (sum(!is.na(fData(a_gse)$GENE)) == length(numerics[numerics==FALSE]))
          successValidation <- TRUE
      }
    }
  }
  return(list("fData"=fData(gse), "successValidation"=successValidation))
}


getDatabaseForMapping <- function(gse) {
  col <- colnames(pData(gse))
  curOrganism <- ""
  for (ch in col) 
    if (grepl("organism", ch)) curOrganism <- tolower(pData(a_gse)[ch][[1]][1])
  
  database <- org.Hs.eg.db
  
  if (curOrganism != "") {
    if (curOrganism == "mus musculus") database <- org.Mm.eg.db
    if (curOrganism == "rattus norvegicus") database <- org.Rn.eg.db
  }
  
  return(database)
}


fitLinearModel <- function(fit, conditions, design, deSize) {
  deList <- list()
  for (i in 1:length(conditions)) {
    
    contrasts <- makeContrasts2(
               c("condition", conditions[[i]]$firstCon, conditions[[i]]$secondCon), 
               levels=design)
    
    fit2 <- contrasts.fit(fit, contrasts)
    fit2 <- eBayes(fit2)
    
    deList[[i]] <- data.table(
                   topTable(fit2, adjust.method="BH", number=deSize, sort.by = "B"), 
                   keep.rownames = T)
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


readDifExprResultsFromFiles <- function() {
  colTypes <- c("character", "character", 
                "double", "double", "double", "double", "double","double")
  files <- dir(path = "./results/dif_expression", 
               full.names = TRUE, recursive = TRUE)
  deList <- list()
  
  for (i in 1:length(files))
    deList[[i]] <- read.table(files[i], header=T, colClasses=colTypes)
  
  return(deList)
}
