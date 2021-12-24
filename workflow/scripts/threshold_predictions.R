suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(plyr)})

main <- function() {
# inputs
  predFile  = {snakemake@input$predictionsSorted}
  thresholdFile = {snakemake@input$thresholds}
  outFile = {snakemake@output$predictionsThresholded}
  method = {snakemake@wildcards$method}
  userThreshold = {snakemake@params$userThresh}
  
  # read in prediction file
  pred = read.table(file=predFile, sep='\t', header=FALSE) %>% 
    setNames(c("chr","start","end","CellType","TargetGene","Score"))
  
  if (userThreshold=="NA" || userThreshold=="None" || is.na(userThreshold)){
    # read in threshold file
    thresholds = read.table(file=thresholdFile, sep='\t', header=TRUE)
    
    if (method=="nearest_tss" || method=="nearest_gene" || method=="within_100kb_of_tss"){
      pred = filter(pred, Score==1)
    } else if (method=="dist_to_tss"){
      t.sub = filter(thresholds, method=="dist_to_tss")
      pred = filter(pred, Score<t.sub$threshold[1])
    } else if (method=="dist_to_gene"){
      t.sub = filter(thresholds, method=="dist_to_gene")
      pred = filter(pred, Score<t.sub$threshold[1])
    } else if (method=="reads_by_dist_to_tss"){
      t.sub = filter(thresholds, method=="reads_by_dist_to_tss")
      pred = filter(pred, Score>t.sub$threshold[1])
    } else if (method=="reads_by_dist_to_tss_norm"){
      t.sub = filter(thresholds, method=="reads_by_dist_to_tss_norm")
      pred = filter(pred, Score>t.sub$threshold[1])
    } else {
      t.sub = filter(thresholds, method==method)
      pred = filter(pred, Score>t.sub$threshold[1])
    }
  } else
  {
    if (method=="dist_to_gene" || method=="dist_to_tss") {
      pred = filter(pred, Score<=userThreshold)
    } else
    {
      pred = filter(pred, Score>=userThreshold)
    }
    
  }
  
  write.table(pred, file=outFile, col.names=FALSE, sep='\t', quote=FALSE, row.names=FALSE)
}

main()