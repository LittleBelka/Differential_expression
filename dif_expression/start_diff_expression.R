library(GEOquery)
# for differentialExpression
source("./dif_expression/couple_universal.R")

options(download.file.method.GEOquery = "libcurl")

dataSetSeries <- c("GSE53986", "GSE50122", "GSE61055", "GSE14308")

for (d in dataSetSeries) {
  message("Start of processing ", d, ".")
  differentialExpression(d)
}