readFile <- function(fileName) {
  con = file(fileName, "r")
  strings <- c()
  i <- 0
  module <- 415
  moduleName <- "GSE81622-GPL10558"
  while (TRUE) {
    line = readLines(con, n = 1)
    if (i != 0) {
      print(i)
      if (length(line) != 0 && grepl(moduleName, line)) {
        line = sub(".*\t", paste(moduleName, '#', module, '\t', sep=''), line)
        module <- module + 1
      }
      strings <- c(strings, line) 
    }
    i <- i + 1
    if (length(line) == 0) break
  }
  close(con)
  return(strings)
}


writeGMT <- function(lines, fileName) {
  if (length(gmt) > 0) {
    fileConn<-file(paste("./some_documents/hs_gmt_files/", fileName, sep=''))
    
    writeLines(lines, fileConn)
    close(fileConn)
  }
}


gmt <- readFile("./some_documents/hs_gmt_files/hs.modules_100001-150000.gmt")
writeGMT(gmt, "hs.modules_10-15.gmt")