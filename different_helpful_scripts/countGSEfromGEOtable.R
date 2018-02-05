library(data.table)

t <- read.table("geo_table.tsv", header=TRUE, sep="\t", fileEncoding="utf-8")

human <- 0
rat <- 0
mouse <- 0
for (i in 1:dim(t)[1]) {
  if (t$Taxonomy[i] == "Homo sapiens") human = human + 1
  if (t$Taxonomy[i] == "Mus musculus") mouse = mouse + 1
  if (t$Taxonomy[i] == "Rattus norvegicus") rat = rat + 1
}

message("Homo sapiens: ", human)
message("Mus musculus: ", mouse)
message("Rattus norvegicus: ", rat)