geneSetEnrichmentAnalysis <- function(deList) {
  gseaStat <- list()
  gseaPlots <- list()
  paths <- list()
  gseaTableResults <- list()
  
  pathways <- getPathways()
  a_pathways <<- pathways

  for (i in 1:length(deList)) {
    ranks <- getRanks(deList[[i]])
    print(head(ranks))
    a_ranks <<- ranks
    
    gseaPlots[[i]] <- plotEnrichment(pathways, ranks, gseaParam = 1)
    gseaStat[[i]] <- calcGseaStat(ranks, na.omit(match(pathways, names(ranks))))
    
    gseaTableResults[[i]] <- getGseaTableResults(pathways, ranks)
    topPathways <- gseaTableResults[[i]][head(order(pval), n=15)][order(NES), pathway]
    plotGseaTable(pathways[topPathways], ranks, gseaTableResults[[i]], gseaParam=0.5)
    
    print(gseaStat[[i]])
    print("_________________________________________________________________-")
  }
  
  resultList <- list("gseaStat"=gseaStat, "gseaPlots"=gseaPlots, 
                     "gseaTableResults"=gseaTableResults)
  return(resultList)
}

getRanks <- function(de) {
  ranksTmp <- de[order(t), list(rn, t)]
  ranks <- as.numeric(unlist(ranksTmp[[2]]))
  ranks <- setNames(ranks, ranksTmp[[1]])
  
  return(ranks)
}

getPathways <- function() {
  mouse.universe <- keys(org.Mm.eg.db, "ENTREZID")
  
  # Selecting reactome gene sets
  pathways <- na.omit(select(reactome.db, keys=mouse.universe, c("PATHID"),
                             keytype = 'ENTREZID'))
  pathways <- split(pathways$ENTREZID, pathways$PATHID)
  
  #pathways <- reactomePathways(names(ranks))
  return(pathways)
}

gseaForDI200 <- function(deList) {
  gseaStat <- list()
  gseaPlots <- list()
  gseaTableResults <- list()
  
  di.file <- "di.dn200.sig.txt"
  diList <- read.table(di.file, header=T)
  
  for (i in 1:length(deList)) {
    deListShort <- getShortListFor200Genes(deList[[i]], diList)
    gseaResultsTmp <- geneSetEnrichmentAnalysisFor200Genes(deListShort)
    gseaStat[[i]] <- gseaResultsTmp$gseaStat 
    gseaPlots[[i]] <- gseaResultsTmp$gseaPlots
    gseaTableResults[[i]] <- gseaTableResults
  }
  
  resultList <- list("gseaStat"=gseaStat, "gseaPlots"=gseaPlots, 
                     "gseaTableResults"=gseaTableResults)
  return(resultList)
}

getShortListFor200Genes <- function(deList, diList) {
  di <- c()
  for (d in diList[[1]]) di <- c(di, d)
  
  deList <- subset(deList, (rn %in% di))
  return(deList)
}

getGseaTableResults <- function(pathways, ranks) {
  fgseaRes <- fgsea(pathways = pathways, stats = ranks,
                    minSize=15,
                    maxSize=500,
                    nperm=1000,
                    nproc=1)
  fgseaRes <- fgseaRes[order(pval)]
  return(fgseaRes)
}

geneSetEnrichmentAnalysisFor200Genes <- function(deList) {
  ranks <- getRanks(deList)
  pathways <- reactomePathways(names(ranks))
  
  a_ranks_200 <<- ranks
  a_path_200 <<- pathways
  
  gseaPlots <- plotEnrichment(pathways, ranks, gseaParam = 1)
  gseaStat <- calcGseaStat(ranks, na.omit(match(pathways, names(ranks))))
  
  gseaTableResults <- getGseaTableResults(pathways, ranks)
  topPathways <- gseaTableResults[head(order(pval), n=15)][order(NES), pathway]
  plotGseaTable(pathways[topPathways], ranks, gseaTableResults, gseaParam=0.5)
  
  print(gseaStat)

  resultList <- list("gseaStat"=gseaStat, "gseaPlots"=gseaPlots, 
                     "gseaTableResults"=gseaTableResults)
  return(resultList)
}

