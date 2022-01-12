suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(plyr)})

main <- function() {
  # snakemake inputs
  predFile  = {snakemake@input$predictionsSorted}
  outFile = {snakemake@output$predictionsThresholded}
  method = {snakemake@wildcards$method}
  threshold = {snakemake@wildcards$threshold} %>% as.numeric()
  
  # read in prediction file
  pred = read.table(file=predFile, sep='\t', header=FALSE) %>% 
    setNames(c("chr","start","end","CellType","TargetGene","Score"))
 
  if (method=="dist_to_tss" || method=="dist_to_gene"){ # methods where high score is worse
    pred.filtered = filter(pred, Score<threshold)
  } else { # all other methods where high score is better
    pred.filtered = filter(pred, Score==threshold || Score>threshold)
  }
  
  write.table(pred.filtered, file=outFile, col.names=FALSE, sep='\t', quote=FALSE, row.names=FALSE)
}

main()
