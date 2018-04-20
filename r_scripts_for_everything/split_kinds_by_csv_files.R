library(dplyr)
library(data.table)

create_csv <- function(newFile, kind) {
  table <- read.csv('some_documents/series_vs_gpl_tsv_with_kind_2.tsv', 
                             sep = '\t', stringsAsFactors = F)
  result <- table %>% 
    filter(GPL != '') %>% 
    filter(Series_taxonomy == kind) %>% 
    select(Series_path, GPL) %>% 
    rename(dataSetSeries = Series_path, availableGPL = GPL)
  
  fwrite(result, newFile, sep = '\t')
  return(result)
}

# create_csv('some_documents/for_slurm_script_2/rat_for_slurm_script_2.tsv', 'rattus norvegicus')
# create_csv('some_documents/for_slurm_script_2/mouse_for_slurm_script_2.tsv', 'mus musculus')
# create_csv('some_documents/for_slurm_script_2/human_for_slurm_script_2.tsv', 'homo sapiens')