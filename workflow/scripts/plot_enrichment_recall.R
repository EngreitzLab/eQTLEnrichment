main <- function() {
  suppressPackageStartupMessages({library(ggplot2)
    library(dplyr)
    library(tidyr)})
  
  # load data
  method = (snakemake@wildcards$method)
  outDir = (snakemake@params$outDir)
  span = (snakemake@params$span) %>% as.character() %>% strsplit(" ") %>% unlist()
  thresholdTableFiles = (snakemake@input$thresholdTables) %>% strsplit(" ") %>% unlist() %>% data.frame() %>% setNames(c('file'))
  outERCurveMaxFull = (snakemake@output$ERCurveMaxFull)
  outERCurveMeanFull = (snakemake@output$ERCurveMeanFull)
  outERCurveMaxZoom = (snakemake@output$ERCurveMaxZoom)
  outERCurveMeanZoom = (snakemake@output$ERCurveMeanZoom)
  outERTable = (snakemake@output$ERCurveTable)
  
  # get list of thresholds
  df = read.table(file=thresholdTableFiles$file[1], header=TRUE, sep="\t") %>% 
    filter(!is.na(prediction.rate.inEnhancer)) %>% unique()
  
  #df = matrix(data=0, nrow=length(span)) %>% as.data.frame()
  
  df$threshold = span %>% as.character()
  df$maxEnrichment = 0
  df$meanEnrichment = 0
  df$recall = 0
  #df = dplyr::select(df, -V1)

  # # format threshold values so that it matches enrichment and count matrix names 
  if (as.numeric(df$threshold[nrow(df)])<10){
    df$threshold[df$threshold=="0"] = "0.0"
    df$threshold[df$threshold=="1"] = "1.0"
    df$threshold[df$threshold=="2"] = "2.0"
    df$threshold[df$threshold=="5"] = "5.0"
    
  }
  df$threshold[df$threshold=="1e+06"] = "1000000"
  df$threshold[df$threshold=="2e+05"] = "200000"
  print(df$threshold)
  
  for (i in 1:nrow(df)){
    # extract threshold from file name of enrichment table 
    threshold.this = df$threshold[i]
    threshold.string = paste0(".", threshold.this, ".tsv")
    enrichmentFile.this = paste0(outDir, "/", method, "/enrichmentTable", threshold.string)
    enrichmentTable.this = read.table(file=enrichmentFile.this, header=TRUE, sep="\t") %>%
      filter(!is.na(enrichment))
    
    # add max and average enrichment or each threshold
    df$maxEnrichment[i] = enrichmentTable.this$enrichment %>% max()
    df$meanEnrichment[i] = enrichmentTable.this$enrichment %>% mean()
  }
  
  # loop through threshold tables for each tissue-biosample pairing
  for (i in 1:nrow(thresholdTableFiles)){
    thresholdTable.this = read.table(file=thresholdTableFiles$file[i], header=TRUE, sep="\t") %>% 
      filter(!is.na(prediction.rate.inEnhancer))
    df$recall = df$recall + thresholdTable.this$prediction.rate.inEnhancer 
  }
  df$recall = df$recall/(nrow(thresholdTableFiles))
  
  # save table
  write.table(df, file=outERTable, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
  
  df$threshold = as.numeric(df$threshold)
  #df$prediction.rate.inEnhancer = df$recall
  
  # graph (color/opacity for threshold?)
  text.mult = 3.5
  
  enr.recall.max = ggplot(df, aes(x=prediction.rate.inEnhancer, y=maxEnrichment, color=threshold)) +
    geom_point() +
    xlab('Recall (fraction of variants in predicted enhancers)') + xlim(c(0,1)) +
    ylab('Maximum enrichment of eQTLs in predicted enhancers') +
    theme_minimal() +
    theme(text=element_text(size=rel(text.mult)), legend.text=element_text(size=rel(2.5)))
  
  enr.recall.mean = ggplot(df, aes(x=prediction.rate.inEnhancer, y=meanEnrichment, color=threshold)) +
    geom_point() +
    xlab('Recall (fraction of variants in predicted enhancers)') + xlim(c(0,1)) +
    ylab('Average enrichment of eQTLs in predicted enhancers') +
    theme_minimal() +
    theme(text=element_text(size=rel(text.mult)), legend.text=element_text(size=rel(2.5)))
  
  enr.recall.max.zoom = ggplot(df, aes(x=prediction.rate.inEnhancer, y=maxEnrichment, color=threshold)) +
    geom_point() +
    xlab('Recall (fraction of variants in predicted enhancers)') + 
    ylab('Maximum enrichment of eQTLs in predicted enhancers') +
    theme_minimal() +
    theme(text=element_text(size=rel(text.mult)), legend.text=element_text(size=rel(2.5)))
  
  enr.recall.mean.zoom = ggplot(df, aes(x=prediction.rate.inEnhancer, y=meanEnrichment, color=threshold)) +
    geom_point() +
    xlab('Recall (fraction of variants in predicted enhancers)') +
    ylab('Average enrichment of eQTLs in predicted enhancers') +
    theme_minimal() +
    theme(text=element_text(size=rel(text.mult)), legend.text=element_text(size=rel(2.5)))
  
  pdf(file=outERCurveMaxFull, width=7, height=5); print(enr.recall.max); dev.off()
  pdf(file=outERCurveMeanFull, width=7, height=5); print(enr.recall.mean); dev.off()
  
  pdf(file=paste0(outERCurveMaxZoom), width=7, height=5); print(enr.recall.max.zoom); dev.off()
  pdf(file=paste0(outERCurveMeanZoom), width=7, height=5); print(enr.recall.mean.zoom); dev.off()
}

main()