library(parallel)
# for differentialExpression
source("./dif_expression/couple_universal.R")
options(download.file.method.GEOquery = "libcurl")

# dataSetSeries <- c("GSE53986", "GSE50122", "GSE61055", "GSE14308")
#dataSetSeries <- c("GSE38230-GPL10123")
dataSetSeries <- c("GSE50122")
fileWithGenes <- "di.dn200.sig.txt"

# difExpr <- function(i) {
#   message("Start of processing ", dataSetSeries[i], ".")
#   differentialExpression(dataSetSeries[i])
# }
# 
# numCores <- detectCores() - 1
# 
# mclapply(1:length(dataSetSeries), difExpr, mc.cores = numCores)


s <- system.time(for (i in 1:length(dataSetSeries)) {
  message("Start of processing ", dataSetSeries[i], ".")
  differentialExpression(dataSetSeries[i], fileWithGenes)
})

# for (i in 1:length(dataSetSeries)) {
#   message("Start of processing ", dataSetSeries[i], ".")
#   differentialExpression(dataSetSeries[i], fileWithGenes)
# }