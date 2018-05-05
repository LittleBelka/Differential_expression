library(data.table)


wrong <- read.table('../some_documents/new_modules_annotations_files/hs.wrong.tsv', 
                    header = T, sep = "\t", colClasses = c("character", "character"))

wrongArr <- c()

for(i in 2:nrow(wrong)) {
  print(i)
  s <- paste(wrong[i,1], '#', wrong[i,2], sep='')
  wrongArr <- c(wrongArr, s)
}


readFile <- function(fileName) {
  con = file(fileName, "r")
  strings <- c()
  j <- 1
  while (TRUE) {
    line = readLines(con, n = 1)
    if (length(line) != 0) {
      module <- gsub("\t.*", '\\1', line)
      if (!(module %in% wrongArr)) {
        strings <- c(strings, line) 
        print(j)
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
fileConn<-file("../some_documents/gmt_files/hs/hs.modules.new.gmt")
writeLines(gmt, fileConn)
close(fileConn)



