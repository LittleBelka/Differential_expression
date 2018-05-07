library(data.table)

annot <- read.csv('./some_documents/new_modules_annotations_files/hs.modules_annotations.tsv', 
                  header = T, sep = "\t")

wrongArr <- c("GSE34200-GPL4091",
              "GSE3892-GPL3255",
              "GSE29023-GPL4091",
              "GSE28329-GPL4091",
              "GSE23631-GPL4091",
              "GSE28331-GPL4091",
              "GSE23633-GPL4091",
              "GSE10878-GPL2879",
              "GSE36825-GPL2879",
              "GSE13239-GPL2879",
              "GSE28114-GPL4091",
              "GSE23005-GPL4091",
              "GSE28329-GPL2873",
              "GSE25893-GPL4091",
              "GSE35735-GPL4091",
              "GSE4003-GPL2913",
              "GSE28331-GPL2873",
              "GSE36823-GPL2879")

annotCorrect <- subset(annot, !(as.character(annot$GSE) %in% wrongArr))

fwrite(annotCorrect, "./some_documents/new_modules_annotations_files/hs.modules_annotations_correct.tsv", sep = '\t')
