library(parallel)
# for differentialExpression
source("./dif_expression/couple_universal.R")
options(download.file.method.GEOquery = "libcurl")

dataSetSeries <- c("GSE53986", "GSE50122", "GSE61055", "GSE14308")
# dataSetSeries <- c("data/GSE50122_series_matrix.txt.gz")
#dataSetSeries <- system('find ../geo/series/ -type f', intern = T)
fileWithGenes <- "di.dn200.sig.txt"

difExpr <- function(i) {
  message("Start of processing ", dataSetSeries[i], ".")
  differentialExpression(dataSetSeries[i], fileWithGenes)
}

numCores <- detectCores() - 1

mclapply(1:length(dataSetSeries), difExpr, mc.cores = numCores)


# for (i in 1:length(dataSetSeries)) {
#   message("Start of processing ", dataSetSeries[i], ".")
#   differentialExpression(dataSetSeries[i], fileWithGenes)
# }