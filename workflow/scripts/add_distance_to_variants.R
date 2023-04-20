suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)})


main <- function() {
  var.file = (snakemake@input$GTExExpressedGenes) 
  TSS.file = (snakemake@params$TSS)
  out.file = (snakemake@output$GTExExpressedGenesWithDistance)
  
  ## merge/join variants with TSS
  # TSS file: no header, chr, start, stop, gene name, score, strand
  TSS = read.table(file=TSS.file, header=FALSE, sep='\t') %>%
    setNames(c("chr", "start", "stop", "gene", "score", "strand"))
  TSS$center = with(TSS, (start + stop)/2)
  
  # variant file: 1-3 (loc), 4 (hgID), 5 (tissue), 6 (gene symbol), 7 (PIP)
  variants = read.table(file=var.file, header=FALSE, sep='\t') %>%
    setNames(c("chr", "start", "stop", "hgID", "tissue", "gene", "PIP", "TPM"))
  
  # add "center" TSS column to variants 
  TSS.col = dplyr::select(TSS, c(gene, center))

  var.merged = inner_join(variants, TSS.col, by="gene")

  ## compute TSS-gene distance
  # subtract and absolute value from variant start loc
  var.merged$distance = with(var.merged, abs(center-start))
  
  var.merged = dplyr::select(var.merged, -center)
  write.table(var.merged, file=out.file, quote=FALSE, sep='\t', col.names=FALSE, row.names=FALSE)
  

  
}


main()
