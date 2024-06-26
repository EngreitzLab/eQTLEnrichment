suppressPackageStartupMessages({library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(data.table)})

# load data
predTableFile = (snakemake@input$predTable)
enrTableFile = (snakemake@input$enrichmentTable)
outERTable = (snakemake@output$ERCurveTable)
method.this = snakemake@wildcards$method

# read in pred table
predTable = fread(predTableFile, header=TRUE, sep="\t")
this.predTable = dplyr::filter(predTable, GTExTissue=="all_matches", Biosample=="all_matches")

# read in enrichment table
enrTable = fread(enrTableFile, header=TRUE, sep="\t")
this.enrTable = dplyr::filter(enrTable, GTExTissue=="all_matches", Biosample=="all_matches") %>% 
	dplyr::select(threshold, enrichment, CI_enr_low, CI_enr_high, SE_log_enr, p_adjust_enr)

# merge by threshold with pred table
df = dplyr::left_join(this.predTable, this.enrTable, by="threshold")
df$method = method.this

# save table
fwrite(df, outERTable, quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)

