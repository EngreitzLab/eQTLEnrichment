suppressPackageStartupMessages({library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)})

# load data
method = (snakemake@wildcards$method)
predTableFile = (snakemake@input$predTable)
enrTableFile = (snakemake@input$enrichmentTable)
outDir = (snakemake@params$outDir)
thresholdSpan = read.table(snakemake@input$thresholdSpan, sep="\t", header=FALSE) %>% setNames("threshold")

GTExTissue.this = (snakemake@wildcards$GTExTissue)
biosample.this = (snakemake@wildcards$biosample)
outERTable = (snakemake@output$ERCurveTable)

# read in prediction table for method x (GTExTissue x biosample)
predTable = read.table(predTableFile, header=TRUE, sep="\t")
predTable$enrichment = 0
predTable$CI_enr_low = 0
predTable$CI_enr_high = 0
predTable$SE_log_enr = 0

# read in enrichment table
enrTable = read.table(enrTableFile, header=TRUE, sep="\t")

# iteratively read in enrichment table per threshold, extract enrichment at that biosample/tissue intersection, add to table
for (i in 1:nrow(thresholdSpan)){
  this.threshold = thresholdSpan$threshold[i]
  this.enrTable = dplyr::filter(enrTable, threshold==this.threshold, GTExTissue==GTExTissue.this, Biosample==biosample.this)
  predTable$enrichment[i] = this.enrTable$enrichment[1]
  predTable$CI_enr_low[i] = this.enrTable$CI_enr_low[1]
  predTable$CI_enr_high[i] = this.enrTable$CI_enr_high[1]
  predTable$SE_log_enr[i] = this.enrTable$SE_log_enr[1]
}

# save table
write.table(predTable, outERTable, quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)

