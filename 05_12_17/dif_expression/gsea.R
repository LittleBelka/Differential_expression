geneSetEnrichmentAnalysis <- function(deList) {
  gseaStat <- list()
  gseaPlots <- list()
  gseaTableResults <- list()
  
  di.file <- "di.dn200.sig.txt"
  diList <- read.table(di.file, header=T)
  diList <- as.character(diList[[1]])
  
  for (i in 1:length(deList)) {
    ranksResults <- getRanks(deList[[i]])
    ranks <- ranksResults$ranks
    
    gseaPlots[[i]] <- plotEnrichment(diList, ranks, gseaParam = 1)
    gseaTableResults[[i]] <- getGseaTableResults(diList, ranks)
    # gseaTableResults[[i]] <- hypergeometricTest(diList, ranksResults$vectorDI)
  }
  
  resultList <- list("gseaPlots"=gseaPlots, 
                     "gseaTableResults"=gseaTableResults)
  return(resultList)
}

getRanks <- function(de) {
  ranksTmp <- de[order(t), list(rn, t)]
  ranks <- as.numeric(unlist(ranksTmp[[2]]))
  ranks <- setNames(ranks, ranksTmp[[1]])
  
  return(list("ranks"=ranks, "vectorDI"=ranksTmp$rn))
}

hypergeometricTest <- function(diList, ranks) {
  overlapCount <- length(intersect(diList, ranks))
  lenList1 <- length(ranks)
  lenList2 <- length(diList)
  
  phyper(overlapCount, lenList1, X-lenList1, lenList2)
  #phyper(184, 20879, 1000, 199) = 0.04005686
  
  return()
}

getGseaTableResults <- function(diList, ranks) {
  fgseaRes <- fgsea(pathways = list(diList), stats = ranks,
                    minSize=15,
                    maxSize=500,
                    nperm=1000,
                    nproc=1)
  fgseaRes <- fgseaRes[order(pval)]
  return(fgseaRes)
}

writeGseaResults <- function(plots, gseaTableResults, conditions, dataSetSeries) {
  dir.create(file.path("./results/", dataSetSeries), showWarnings = FALSE)
  filePath <- paste("./results/", dataSetSeries, "/", sep="")
  
  dir.create(file.path(filePath, "fgsea"), showWarnings = FALSE)
  dir.create(file.path(filePath, "plots"), showWarnings = FALSE)
  
  for (i in 1:length(conditions)) {
    nameFileWithTableResults <- paste(filePath,
                              "/fgsea/",        
                              conditions[[i]]$firstCon, ".vs.", 
                              conditions[[i]]$secondCon, ".tsv", sep="")
    
    nameFileWithPlot <- paste(filePath,
                      "/plots/",        
                      conditions[[i]]$firstCon, ".vs.", 
                      conditions[[i]]$secondCon, ".png", sep="")
    
    fwrite(gseaTableResults[[i]], nameFileWithTableResults)
    ggsave(filename=nameFileWithPlot, plot=plots[[i]], 
           width = 9, height = 6, units = "cm")
  }
}
