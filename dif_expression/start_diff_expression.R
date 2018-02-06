library(parallel)
# for differentialExpression
source("./dif_expression/couple_universal.R")
# options(download.file.method.GEOquery = "libcurl")

downloadLibrariesAndSources()

# gseVSgpl <- findGSEvsGPL()
# print(gseVSgpl)

# gseVSgpl <- data.frame("dataSetSeries"="data_3/GSE57802_series_matrix.txt.gz",
#                  "availableGPL"="data_3/GPL13158.tsv", stringsAsFactors = F)

# gseVSgpl <- data.frame("dataSetSeries"="data_4/GSE57378_series_matrix.txt.gz",
#                        "availableGPL"="data_4/GPL10558.tsv", stringsAsFactors = F)

# gseVSgpl <- data.frame("dataSetSeries"="data/GSE61055_series_matrix.txt.gz",
#                        "availableGPL"="data/GPL6887.tsv", stringsAsFactors = F)

# gseVSgpl <- data.frame("dataSetSeries"="data_4/GSE38230-GPL10558_series_matrix.txt.gz",
#                        "availableGPL"="data_4/GPL10558.tsv", stringsAsFactors = F)

# gseVSgpl <- data.frame("dataSetSeries"="data_4/GSE49968_series_matrix.txt.gz",
#                        "availableGPL"="data_4/GPL7202.tsv", stringsAsFactors = F)

# gseVSgpl <- data.frame("dataSetSeries"="data_4/GSE7033_series_matrix.txt.gz",
#                        "availableGPL"="data_4/GPL3050.tsv", stringsAsFactors = F)

# gseVSgpl <- data.frame("dataSetSeries"="data_4/GSE29262_series_matrix.txt.gz",
#                        "availableGPL"="data_4/GPL1261.tsv", stringsAsFactors = F)

# gseVSgpl <- data.frame("dataSetSeries"="data_4/GSE35088_series_matrix.txt.gz",
#                        "availableGPL"="data_4/GPL6101.tsv", stringsAsFactors = F)

gseVSgpl <- data.frame("dataSetSeries"="data_4/GSE28026_series_matrix.txt.gz",
                       "availableGPL"="data_4/GPL570.tsv", stringsAsFactors = F)

# series <- c("data_4/GSE49968_series_matrix.txt.gz", "data_4/GSE7033_series_matrix.txt.gz", 
#             "data_4/GSE29262_series_matrix.txt.gz")
# 
# platforms <- c("data_4/GPL7202.tsv", "data_4/GPL3050.tsv", "data_4/GPL1261.tsv")
# 
# gseVSgpl <- data.frame("dataSetSeries"=series, "availableGPL"=platforms, stringsAsFactors = F)

#dataSetSeries <- c("GSE53986", "GSE50122", "GSE61055", "GSE14308")
# dataSetSeries <- c("GSE10759")

fileWithGenes <- "some_documents/di.dn200.sig.txt"

difExpr <- function(i) {
  message("Start of processing ", gseVSgpl$dataSetSeries[i], ".")
  differentialExpression(
        gseVSgpl$dataSetSeries[i], gseVSgpl$availableGPL[i], fileWithGenes)
}

numCores <- detectCores() - 1

mclapply(1:nrow(gseVSgpl), difExpr, mc.cores = numCores)


# for (i in 1:length(dataSetSeries)) {
#   message("Start of processing ", dataSetSeries[i], ".")
#   differentialExpression(dataSetSeries[i], fileWithGenes)
# }
