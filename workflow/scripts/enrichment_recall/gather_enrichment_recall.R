suppressPackageStartupMessages({library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)})

# load data
predTableFile = (snakemake@input$predTable)
enrTableFile = (snakemake@input$enrichmentTable)

GTExTissue.this = (snakemake@wildcards$GTExTissue)
biosample.this = (snakemake@wildcards$biosample)
outERTable = (snakemake@output$ERCurveTable)

# read in prediction table for method x (GTExTissue x biosample)
predTable = read.table(predTableFile, header=TRUE, sep="\t")

# read in enrichment table
enrTable = read.table(enrTableFile, header=TRUE, sep="\t")
enrTable = read.table(enrTableFile, header=TRUE, sep="\t")
this.enrTable = dplyr::filter(enrTable, GTExTissue==GTExTissue.this, Biosample==biosample.this) %>% 
	dplyr::select(threshold, enrichment, CI_enr_low, CI_enr_high, SE_log_enr)

# merge by threshold with pred table
predTable = dplyr::left_join(predTable, this.enrTable, by="threshold")

# save table
write.table(predTable, outERTable, quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)

