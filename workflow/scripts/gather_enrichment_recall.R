suppressPackageStartupMessages({library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)})

# load data
method = (snakemake@wildcards$method)
method.slash = paste0("/", method, "/")
predTableFile = (snakemake@input$predTable)
enrFiles = (snakemake@input$enrichmentMatrices)%>% strsplit(" ") %>% unlist() %>% data.frame() %>% setNames("file")
print(enrFiles)
enrFiles = dplyr::filter(enrFiles, grepl(method.slash, file, fixed=TRUE))
print(enrFiles)
outDir = (snakemake@params$outDir)
thresholdSpan = (snakemake@params$thresholdSpan) %>% as.character() %>% strsplit(" ") %>% unlist() %>% as.numeric()
GTExTissue.this = (snakemake@wildcards$GTExTissue)
biosample.this = (snakemake@wildcards$Biosample)
outERTable = (snakemake@output$ERCurveTable)

# read in prediction table for method x (GTExTissue x biosample)
predTable = read.table(predTableFile, header=TRUE, sep="\t")
predTable$enrichment = 0

# iteratively read in enrichment table per threshold, extract enrichment at that biosample/tissue intersection, add to table
for (i in 1:length(thresholdSpan)){
  #this.threshold = thresholdSpan[i]
  #this.fileName = paste0("enrichmentTable.thresh", this.threshold, ".tsv")
  #this.enrFile = file.path(outDir, method, "enrichmentTables", this.fileName)
  this.enrFile = enrFiles$file[i]
  this.enrTable = read.table(this.enrFile, header=TRUE, sep="\t")
  this.enrTable = dplyr::filter(this.enrTable, GTExTissue==GTExTissue.this, Biosample==biosample.this)
  predTable$enrichment[i] = this.enrTable$enrichment[1]
}

# save table
write.table(predTable, outERTable, quote=FALSE, sep="\t", col.names=TRUE)

