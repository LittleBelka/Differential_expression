library(stringr)

con = file('some_documents/modules_annot_txt_files/hs.modules_annotations.txt', "r")
strings <- c()
i <- 1
while (TRUE) {
  line = readLines(con, n = 1)
  if (length(line) != 0) {
    if (grepl('[^-]GPL', line)) {
      line <- str_replace(line, '(GPL\\d+)\\1(\\d+)', '\\1#\\2\t')
      message(i, ' ', line)
    }
    strings <- c(strings, line) 
    i <- i + 1
  }
  if (length(line) == 0) break
}
close(con)


# fileConn<-file('some_documents/modules_annot_txt_files/hs.modules_annotations_2.txt')
# writeLines(strings, fileConn)
# close(fileConn)

