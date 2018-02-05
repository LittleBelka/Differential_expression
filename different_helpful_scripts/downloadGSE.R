library(GEOquery)
library(data.table)

t <- read.table("geo_table.tsv", header=TRUE, sep="\t", fileEncoding="utf-8")

human <- 0
rat <- 0
mouse <- 0
currentHuman <- ""
currentMouse <- ""
currentRat <- ""
existedHuman <- array(list.files("./data1/human/downloaded_matrices", pattern = "*.gz"))
existedMouse <- array(list.files("./data1/mouse/downloaded_matrices", pattern = "*.gz"))
existedRat <- array(list.files("./data1/rat/downloaded_matrices", pattern = "*.gz"))

# for (i in 1:dim(t)[1]) {
for (i in 1:5) {
  if (t$Taxonomy[i] == "Homo sapiens") {
    currentHuman <- as.character(t$Accession[i])
    z <- sapply(existedHuman, function(x) grepl(currentHuman, x))
    
    if (length(z[z==TRUE]) == 0) {
      tryCatch(getGEO(currentHuman, destdir = "./data2/human/downloaded_matrices")[[1]],
               error = function(e) message(e))
      human = human + 1
    }
  }
  
  if (t$Taxonomy[i] == "Mus musculus") {
    currentMouse <- as.character(t$Accession[i])
    z <- sapply(existedMouse, function(x) grepl(currentMouse, x))
    
    if (length(z[z==TRUE]) == 0) {
      tryCatch(getGEO(currentMouse, destdir = "./data2/mouse/downloaded_matrices")[[1]],
               error = function(e) message(e))
      mouse = mouse + 1
    }
  }
  
  if (t$Taxonomy[i] == "Rattus norvegicus") {
    currentRat <- as.character(t$Accession[i])
    z <- sapply(existedRat, function(x) grepl(currentRat, x))
    
    if (length(z[z==TRUE]) == 0) {
      tryCatch(getGEO(currentRat, destdir = "./data2/rat/downloaded_matrices")[[1]],
               error = function(e) message(e))
      rat = rat + 1
    }
  }
}

message("Homo sapiens: ", human, ", current GSE: ", currentHuman)
message("Mus musculus: ", mouse, ", current GSE: ", currentMouse)
message("Rattus norvegicus: ", rat, ", current GSE: ", currentRat)
