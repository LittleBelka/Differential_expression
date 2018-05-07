readFile <- function(fileName) {
  con = file(fileName, "r")
  strings <- c()
  j <- 1
  k <- 1
  while (TRUE) {
    line = readLines(con, n = 1)
    if (length(line) != 0) {
      module <- gsub("\t.*", '\\1', line)
      genes <- sub(".*\t", '', line)
      count <- strsplit(genes,split=',', fixed=TRUE)
      if (length(count[[1]]) == 7000 || length(count[[1]]) == 200) {
        print(k)
        k <- k + 1
      } else {
        t <- paste(module, '\t', length(count[[1]]), sep='')
        strings <- c(strings, t) 
        message(j, ' ', module, ' ', length(count[[1]]))
        j <- j + 1
      }
    }
    if (length(line) == 0) break
  }
  close(con)
  return(strings)
}

print("Read gmt file")
gmt <- readFile('../some_documents/gmt_files/hs/hs.modules.gmt')

print("Write new gmt")
fileConn<-file("../some_documents/gmt_files/hs/hs.modules_wrong.txt")
writeLines(gmt, fileConn)
close(fileConn)