library(GEOquery)
library(limma)
library(org.Mm.eg.db)
library(gtools)
library(stringr)
library(fgsea)
library(data.table)
library(reactome.db)
# for collapseBy
source("https://raw.githubusercontent.com/assaron/r-utils/master/R/exprs.R")
# for conditions handling
source("./difExpression/couple_universal_handler.R")
# for gene set enrichment analysis
source("./difExpression/gsea.R")

options(download.file.method.GEOquery = "libcurl")
#gse <- getGEO("GSE53986", destdir = ".")[[1]]
#gse <- getGEO("GSE50122", destdir = ".")[[1]] 
gse <- getGEO("GSE61055", destdir = ".")[[1]] 
#gse <- getGEO("GSE14308", destdir = ".")[[1]] 

characteristics <- getCharacteristicsColumns(gse)
conditionLists <- list()

if (length(characteristics) != 0) {
  message("There are characteristics columns in the samples table.")
  print(characteristics)
  conditionLists <- getConditionsFromCharacteristics(gse, characteristics)
}
if (length(conditionLists) > 0) {
  pData(gse)$condition <- fillGseConditionColumn(conditionLists)
  message("Column 'condition' was created and filled in the samples table.")
} else {
  message("Characteristics columns in the samples table don't exist or were unhelpful.
          So start to parse title column.")
  getConditionsFromTitle()
}

con <- pData(gse)$condition

fData(gse) <- provideValidOfSomeColumns(gse)

es <- collapseBy(gse, fData(gse)$ENTREZ_GENE_ID, FUN=median)
es <- es[!grepl("///", rownames(es)), ]
es <- es[rownames(es) != "", ]

message("Garbage was deleted from gene table.")

fData(es) <- data.frame(row.names = rownames(es))

fData(es)$symbol <- mapIds(org.Mm.eg.db, keys = rownames(es), column = "SYMBOL", keytype = "ENTREZID")

exprs(es) <- normalizeBetweenArrays(log2(exprs(es)+1), method="quantile")

es.design <- model.matrix(~0+condition, data=pData(es))

fit <- lmFit(es, es.design)
conditions <- getConditionsForBuildingLinearModel(pData(gse)$condition)

message("Conditions combinations were received for filling of contrast matrix.")

# Show received pairs of comparisons
for (i in 1:length(conditions)) 
  message(conditions[[i]][1], " ", conditions[[i]][2])

deSize <- dim(fData(es))[1]
deList <- fitLinearModel(fit, conditions, es.design, deSize)
message("Linear Models were fitted and saved in 'deList'.")

writeDifExprResultsToFiles(deList, conditions)
message("Linear Models were written to files.")

#deList <- readDifExprResultsFromFiles()
#message("Linear Models were read from files and stored in 'deList'.")

gseaResults <- geneSetEnrichmentAnalysis(deList)
gseaStat <- gseaResults$gseaStat
plots <- gseaResults$gseaPlots
gseaTableResults <- gseaResults$gseaTableResults
message("Gene set enrichment analysis was done.")
