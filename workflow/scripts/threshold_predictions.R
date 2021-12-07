suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(plyr)})

main <- function() {
# inputs
  predFile  = {snakemake@input$predictionsSorted}
  outFile = {snakemake@output$predictionsThresholded}
  method = {snakemake@wildcards$method}
  
# read in prediction file
  pred = read.table(file=predFile, sep='\t', header=FALSE) %>% 
    setNames(c("chr","start","end","CellType","TargetGene","Score"))
  
if (method=="nearest_tss" || method=="nearest_gene" || method=="within_100kb_of_tss"){
  pred = filter(pred, Score==1)
}

if (method=="dist_to_tss"){
  pred = filter(pred, Score<25000)
}
  
if (method=="dist_to_gene"){
  pred = filter(pred, Score<54000)
}
  
if (method=="reads_by_dist_to_tss"){
  pred = filter(pred, Score>0.00625)
}
  
if (method=="reads_by_dist_to_tss_norm"){
  pred = filter(pred, Score>0.00625)
  
}
  
  write.table(pred, file=outFile, col.names=FALSE, sep='\t', quote=FALSE, row.names=FALSE)
}

main()
