library(data.table)


annot <- read.table('./some_documents/hs_gmt_files/hs.modules_annotations_300001-317340.tsv', 
                       header = T, sep = "\t", colClasses = c("character", "character","character", "character"))

module <- 253
moduleName <- "GSE78061-GPL6244"

for(i in 2:nrow(annot)) {
  if (moduleName == as.character(annot[i,1])) {
    print(i)
    annot[i,3] = module
    module = module + 1
  }
}

fwrite(annot[2:nrow(annot),], "./some_documents/hs_gmt_files/annot_30-317.tsv", sep = '\t', col.names=F)

