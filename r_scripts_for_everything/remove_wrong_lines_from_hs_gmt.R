library(data.table)


wrong <- read.table('some_documents/new_modules_annotations_files_2/mm.wrong.tsv', 
                    header = T, sep = "\t", colClasses = c("character", "character"))

wrongArr <- c()

for(i in 2:nrow(wrong)) {
  print(i)
  s <- paste(wrong[i,1], '#', wrong[i,2], sep='')
  wrongArr <- c(wrongArr, s)
}

# wrongArr <- c("GSE34200-GPL4091",
#               "GSE3892-GPL3255",
#               "GSE29023-GPL4091",
#               "GSE28329-GPL4091",
#               "GSE23631-GPL4091",
#               "GSE28331-GPL4091",
#               "GSE23633-GPL4091",
#               "GSE10878-GPL2879",
#               "GSE36825-GPL2879",
#               "GSE13239-GPL2879",
#               "GSE28114-GPL4091",
#               "GSE23005-GPL4091",
#               "GSE28329-GPL2873",
#               "GSE25893-GPL4091",
#               "GSE35735-GPL4091",
#               "GSE4003-GPL2913",
#               "GSE28331-GPL2873",
#               "GSE36823-GPL2879")

readFile <- function(fileName) {
  con = file(fileName, "r")
  strings <- c()
  j <- 1
  while (TRUE) {
    line = readLines(con, n = 1)
    if (length(line) != 0) {
      module <- gsub("\t.*", '\\1', line)
      # module <- sub("#.*", '', module)
      # module <- sub('_', '-', module)
      # strings <- c(strings, module) 
      # module <- gsub("#.*\t.*", '\\1', line)
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
gmt <- readFile('some_documents/gmt_files/mm.modules.2.gmt')

# gmt_genequery <- readFile('some_documents/genequery_gmt/mm.modules.gmt')
# 
# gmt_rt <- readFile('some_documents/gmt_files/rt.modules.gmt')
# gmt_genequery_rt <- readFile('some_documents/genequery_gmt/rt.modules.gmt')

print("Write new gmt")
fileConn<-file("some_documents/gmt_files/mm.modules.new.gmt")
writeLines(gmt, fileConn)
close(fileConn)



