suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(data.table)})


main <- function() {
  biosamples = (snakemake@params$biosamples) %>% as.character() %>% strsplit(" ") %>% unlist() 
  tissues = (snakemake@params$tissues) %>% as.character() %>% strsplit(" ") %>% unlist() 
  out.file = (snakemake@output$map)

  df = data.frame(biosample = biosamples, tissue = tissues)

  write.table(df, out.file, quote=FALSE, sep='\t', col.names=TRUE, row.names=FALSE)
  
}


main()
