suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)})


main <- function() {
  var.file = (snakemake@input$filteredGTExVariants) 
  distances = (snakemake@params$distances) %>% unlist() %>% as.numeric()
  TSS.file = (snakemake@params$TSS)
  out.files = (snakemake@output$filteredGTExVariantsOut) %>% strsplit(" ") %>% unlist()
  
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

  ## filter for each distance and write output file without distance col
  # for first threshold
  filtered.this = dplyr::filter(var.merged, var.merged$distance<=distances[1])
  out.this = dplyr::select(filtered.this, -c(distance, center))
  write.table(out.this, file=out.files[1], quote=FALSE, sep='\t', col.names=FALSE, row.names=FALSE)
  
  for (i in 2:length(distances)) {
    filtered.this = dplyr::filter(var.merged, var.merged$distance<=distances[i], var.merged$distance>=distances[i-1])
    out.this = dplyr::select(filtered.this, -c(distance, center))
    write.table(out.this, file=out.files[i], quote=FALSE, sep='\t', col.names=FALSE, row.names=FALSE)
  }

  
}


main()
