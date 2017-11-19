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
  
  for (ch in characteristics) {
    uniqChar <- unique(pData(gse)[ch])
    
    if (length(rownames(uniqChar)) > 1) 
      conditionsList <- getConditionsFromUniqValues(
                        gse, pData(gse)[ch], conditionsList, uniqChar)
  }
  return(conditionsList)
}

getConditionsFromUniqValues <- function(gse, charColumns, conditionsList, uniqChar) {
  values <- sub("^.*: ", "", charColumns[[1]])
  numerics <- is.na(suppressWarnings(as.numeric(values)))
  
  if (length(numerics[numerics==FALSE]) == 0) {
    tmpCondition <- c()
    complicatedCondition <- FALSE
    tmpComplicatedCondition <- createTmpComplicatedConditionList(uniqChar)
    
    for (v in values) {
      if (grepl("\\(.*\\)", v)) 
        v <- str_replace(v, '.*\\((.*)\\).*', '\\1')
      if (isTooComplicatedCondition(v)) complicatedCondition <- TRUE
      
      if (complicatedCondition) {
        resultCondition <- handleComplicatedCondition(tmpComplicatedCondition, v)
        v <- resultCondition$condition
        tmpComplicatedCondition <- resultCondition$tmpComplicatedCondition
      } else 
        v <- str_replace_all(v, "[\\+ _\\-:;\\.\\*!'\\/]*", "")
      
      tmpCondition <- c(tmpCondition, v)
    }
    conditionsList[[length(conditionsList)+1]] <- tmpCondition
  }
  return(conditionsList)
}

createTmpComplicatedConditionList <- function(uniqChar) {
  tmpComplicatedCondition <- list()
  for (i in 1:length(rownames(uniqChar))) {
    tmpComplicatedCondition[[i]] <- list("condition"="", 
                                         replacement=intToUtf8(64+i))
  }
  return(tmpComplicatedCondition)
}

handleComplicatedCondition <- function(tmpComplicatedCondition, v) {
  newConditionName <- isContainedInConditionList(tmpComplicatedCondition, v)
  if (newConditionName != "") {
    v <- newConditionName 
  } else {
    newConditionIndex <- getFreeListPosition(tmpComplicatedCondition)
    tmpComplicatedCondition[[newConditionIndex]]$condition = v
    v = tmpComplicatedCondition[[newConditionIndex]]$replacement
  }
  return(list("condition"=v, "tmpComplicatedCondition"=tmpComplicatedCondition))
}

getFreeListPosition <- function(tmpComplicatedCondition) {
  i = 0
  found <- FALSE
  while(!found && length(tmpComplicatedCondition) > i) {
    i = i + 1
    if (tmpComplicatedCondition[[i]]$condition == "") 
      found <- TRUE
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
  if (!("ENTREZ_GENE_ID" %in% names(fData(gse)))) {
    columnNames <- names(fData(gse))
    
    for (i in 1:length(columnNames)) {
      if (sapply(columnNames[i], tolower)[[1]] == "entrez_gene_id")
        names(fData(gse))[names(fData(gse)) == columnNames[i]] <- "ENTREZ_GENE_ID"
    }
  }
  return(fData(gse))
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

writeDifExprResultsToFiles <- function(deList, conditions) {
  for (i in 1:length(conditions)) {
    nameFile <- paste("./results/text_data/difExpr/",
                      conditions[[i]]$firstCon, ".vs.", 
                      conditions[[i]]$secondCon, ".tsv", sep="")
    write.table(deList[[i]], nameFile, sep="\t", quote = F, row.names = F)
  }
}

readDifExprResultsFromFiles <- function() {
  colTypes <- c("character", "character", 
                "double", "double", "double", "double", "double","double")
  files <- dir(path = "./results/text_data/difExpr", 
               full.names = TRUE, recursive = TRUE)
  deList <- list()
  
  for (i in 1:length(files))
    deList[[i]] <- read.table(files[i], header=T, colClasses=colTypes)
  
  return(deList)
}
