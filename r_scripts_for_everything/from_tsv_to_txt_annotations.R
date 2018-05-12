library(data.table)
library(dplyr)
library(stringr)


writeToFile <- function(t, filename) {
  fileConn<-file(filename)
  toWrite <- c()
  
  for (i in 1:length(t$Module)) {
    s <- paste(t$Module[[i]], t$Module_name[[i]], sep='\t')
    toWrite <- c(toWrite, s)
  }
  
  writeLines(toWrite, fileConn)
  close(fileConn)
}


# gseTable <- read.csv('some_documents/new_modules_annotations_files/rt.modules_annotations.tsv',
#                        header = T, sep = "\t", fill = TRUE)

# gseTable <- read.csv('some_documents/new_modules_annotations_files/mm.modules_annotations.tsv',
#                        header = T, sep = "\t", fill = TRUE)

gseTable <- read.csv('some_documents/new_modules_annotations_files/hs.modules_annotations.tsv',
                       header = T, sep = "\t", fill = TRUE)

result <- gseTable %>% 
  mutate(Module_name = str_replace_all(Module_name, '\t', '')) %>%
  mutate(Module = paste(GSE, '#', Module_number, sep = ''))  %>% 
  select(Module, Module_name)


# writeToFile(result, 'some_documents/modules_annot_txt_files/rt.modules_annotations.txt')
# writeToFile(result, 'some_documents/modules_annot_txt_files/mm.modules_annotations.txt')
writeToFile(result, 'some_documents/modules_annot_txt_files/hs.modules_annotations.txt')