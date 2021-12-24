suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(plyr)
  library(stringr)})

main <- function() {
# inputs
  thresholdFiles  = {snakemake@input$thresholdTables} %>% as.character() %>% strsplit(", ") %>% unlist() %>% data.frame() %>% setNames(c("fileName"))
  outFile = {snakemake@output$thresholds}
  rate = {snakemake@params$overlapRate}
  methods = {snakemake@params$methods}

  # iterate through the four methods
  thresholds = data.frame(matrix(data=0, nrow=length(methods), ncol=2)) %>% setNames(c("method", "threshold"))
  thresholds$method = methods

  for (i in 1:nrow(thresholds)){
    # filter list of files to this method
    tables.this = filter(thresholdFiles, str_detect(string=fileName, pattern=paste0("/",thresholds$method[i], "/")))
    # loop through files (each a different tissue)
    for (j in 1:nrow(tables.this)){
      df = read.table(tables.this$fileName[j], header=TRUE, sep="\t")
      # find cut-off value for [rate]
      df = df[order(df$prediction.rate.inEnhancer),] # sort df in ascending order of prediction.rate.inEnhancer
      df.filtered = filter(df, prediction.rate.inEnhancer>=rate)  # filter to rate above desired value
      #if there are no thresholds with a rate above the value
      # choose the threshold with the highest rate (last in the df)
      # otherwise, choose the lowest threshold in the filtered df
      threshold.this = dplyr::if_else(nrow(df.filtered)==0,
                                    df$threshold[nrow(df)], 
                                    df.filtered$threshold[1])
      # add to average
      thresholds$threshold[i] = thresholds$threshold[i] + threshold.this
    }
    thresholds$threshold[i] = thresholds$threshold[i]/nrow(tables.this)
  }
  
  write.table(thresholds, file=outFile, col.names=TRUE, sep='\t', quote=FALSE, row.names=FALSE)
}

main()
