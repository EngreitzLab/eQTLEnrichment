main <- function() {
  suppressPackageStartupMessages({library(ggplot2)
    library(dplyr)
    library(tidyr)})
  
  # load data
  enrichmentTableFiles = (snakemake@input$enrichmentTables) %>% strsplit(" ") %>% unlist() %>% data.frame() %>% setNames(c('file'))
  method = (snakemake@wildcards$method)
  thresholdTableFiles = (snakemake@input$thresholdTables) %>% strsplit(" ") %>% unlist() %>% data.frame() %>% setNames(c('file'))
  outERCurveMaxFull = (snakemake@output$ERCurveMaxFull)
  outERCurveMeanFull = (snakemake@output$ERCurveMeanFull)
  outERCurveMaxZoom = (snakemake@output$ERCurveMaxZoom)
  outERCurveMeanZoom = (snakemake@output$ERCurveMeanZoom)
  outERTable = (snakemake@output$ERCurveTable)
  
  # get list of thresholds
  df = read.table(file=thresholdTableFiles$file[1], header=TRUE, sep="\t") %>% 
    filter(!is.na(prediction.rate.inEnhancer)) %>% unique()
  df$maxEnrichment = 0
  df$meanEnrichment = 0
  df$recall = 0
  
  # # format threshold values so that it matches enrichment and count matrix names 
  # if (method=='dist_to_tss' || method=='dist_to_gene'){
  #   df$threshold = format(df$threshold, scientific=FALSE, trim=TRUE) %>% as.character()
  #   df$threshold[df$threshold=="100000"] = "100000.0"
  # } else if (method=='reads_by_dist_to_tss' || method=='reads_by_dist_to_tss_norm'){
  #   df$threshold= format(df$threshold, scientific=FALSE) %>% as.character()
  #   df$threshold = sub("0+$", "", df$threshold) 
  # }
  
  for (i in 1:nrow(df)){
    # extract threshold from file name of enrichment table 
    threshold.this = df$threshold[i]
    threshold.string = paste0(".", threshold.this, ".tsv")
    enrichmentFile.this = filter(enrichmentTableFiles, grepl(threshold.string, enrichmentTableFiles$file))
    enrichmentTable.this = read.table(file=enrichmentFile.this$file[1], header=TRUE, sep="\t") %>%
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