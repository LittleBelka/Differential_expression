create_csv <- function(newFile, kind) {
  table <- read.csv('some_documents/series_vs_gpl_tsv_with_kind.csv', 
                             sep = '\t', stringsAsFactors = F)
  result <- table %>% 
    filter(GPL != '') %>% 
    filter(Series_taxonomy == kind) %>% 
    select(Series_path, GPL) %>% 
    rename(dataSetSeries = Series_path, availableGPL = GPL)
  
  
  fwrite(result, newFile, sep = '\t')
}

# create_csv('some_documents/rat_for_slurm_script.csv', 'rattus norvegicus')
create_csv('some_documents/mouse_for_slurm_script.csv', 'mus musculus')
create_csv('some_documents/human_for_slurm_script.csv', 'homo sapiens')