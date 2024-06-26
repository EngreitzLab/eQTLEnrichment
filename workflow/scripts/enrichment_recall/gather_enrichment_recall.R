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
this.predTable = dplyr::filter(predTable, GTExTissue==GTExTissue.this, Biosample==biosample.this)

# read in enrichment table
enrTable = read.table(enrTableFile, header=TRUE, sep="\t")
this.enrTable = dplyr::filter(enrTable, GTExTissue==GTExTissue.this, Biosample==biosample.this) %>% 
	dplyr::select(GTExTissue, Biosample, threshold, enrichment, CI_enr_low, CI_enr_high, SE_log_enr, p_adjust_enr)

# merge by threshold with pred table
df = dplyr::left_join(this.predTable, this.enrTable, by=c("threshold", "GTExTissue", "Biosample"))

# save table
write.table(df, outERTable, quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)

