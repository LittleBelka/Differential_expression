library(parallel)
# for differentialExpression
source("./dif_expression/couple_universal.R")
# options(download.file.method.GEOquery = "libcurl")


# args <- commandArgs(trailingOnly = TRUE)
# dataSetSeries <- args[1]
# gpl <- args[2]
# 
downloadLibrariesAndSources()
# 
# message(dataSetSeries, ' ', gpl)
# differentialExpression(dataSetSeries, gpl)


# gseVSgpl <- findGSEvsGPL()
# print(gseVSgpl)

gseVSgpl <- data.frame("dataSetSeries"="data_3/GSE57802_series_matrix.txt.gz",
                       "availableGPL"="data_3/GPL13158.tsv", stringsAsFactors = F)

# gseVSgpl <- data.frame("dataSetSeries"="data/GSE61055_series_matrix.txt.gz",
#                        "availableGPL"="data/GPL6887.tsv", stringsAsFactors = F)

# gseVSgpl <- data.frame("dataSetSeries"="data_4/GSE37815_series_matrix.txt.gz",
#                        "availableGPL"="data_4/GPL6102.tsv", stringsAsFactors = F)
                       

# gseVSgpl <- data.frame("dataSetSeries"="data_4/GSE28026_series_matrix.txt.gz",
#                        "availableGPL"="data_4/GPL570.tsv", stringsAsFactors = F)

fileWithGenes <- "some_documents/di.dn200.sig.txt"

difExpr <- function(i) {
 message("Start of processing ", gseVSgpl$dataSetSeries[i], ".")
 differentialExpression(
       gseVSgpl$dataSetSeries[i], gseVSgpl$availableGPL[i])
}

numCores <- detectCores() - 1

mclapply(1:nrow(gseVSgpl), difExpr, mc.cores = numCores)


# for (i in 1:length(dataSetSeries)) {
#   message("Start of processing ", dataSetSeries[i], ".")
#   differentialExpression(dataSetSeries[i], fileWithGenes)
# }
