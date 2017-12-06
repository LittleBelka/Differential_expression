library(parallel)
# for differentialExpression
source("./dif_expression/couple_universal.R")
options(download.file.method.GEOquery = "libcurl")

dataSetSeries <- c("GSE53986", "GSE50122", "GSE61055", "GSE14308")

difExpr <- function(i) {
  message("Start of processing ", dataSetSeries[i], ".")
  differentialExpression(dataSetSeries[i])
}

numCores <- detectCores() - 1

mclapply(1:length(dataSetSeries), difExpr, mc.cores = numCores)


