suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(plyr)})

main <- function() {
# inputs
  predFile  = {snakemake@input$predictionsSorted}
  outFile = {snakemake@input$predictionsThresholded}
  method = {snakemake@wildcards$method}
  
# read in prediction file
  pred = read.table(file=predFile, sep='\t', header=FALSE) %>% 
    setNames(c("chr","start","end","CellType","TargetGene","Score"))
  
if (method=="nearest_tss" || method=="nearest_gene" || method=="within_100kb_of_tss"){
  pred = filter(pred, Score==1)
}

if (method=="dist_to_tss"){
  
}
  
if (method=="dist_to_gene"){
    
}
  
if (method=="reads_by_dist_to_tss"){
  
}
  
if (method=="reads_by_dist_to_tss_norm"){
  
}
  
  write.table(pred, file=outFile, header=FALSE, sep='\t', quote=FALSE, row.names=FALSE)
}

main()
