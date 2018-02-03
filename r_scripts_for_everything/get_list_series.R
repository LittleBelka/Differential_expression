dataSetSeries <- system('find ../geo/series/ -type f', intern = T)

fileConn <- file("series.txt")
writeLines(dataSetSeries, fileConn)
close(fileConn)

listGSE <- c()
for (i in 1:length(dataSetSeries)) {
  tmp <- sub('(.*/)*', '', dataSetSeries[i])
  tmp <- sub('_series_matrix.txt.gz', '', tmp)
  if (tmp == 'GSE61055' || tmp == 'GSE53986' 
      || tmp == 'GSE14308' || tmp == 'GSE50122'
      || tmp == 'GSE9943' || tmp == 'GSE38230') { 
    listGSE <- c(listGSE, dataSetSeries[i])
    print(dataSetSeries[i])
  }
}

fileConn <- file("list_series_for_checking.txt")
writeLines(listGSE, fileConn)
close(fileConn)