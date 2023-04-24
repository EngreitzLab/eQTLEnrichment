suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)})


main <- function() {
  var.file = (snakemake@input$GTExExpressedGenes) 
  TSS.file = (snakemake@params$TSS)
  out.file = (snakemake@output$GTExExpressedGenesWithDistance)
  distances = (snakemake@params$distances) %>% as.character() %>% strsplit(" ") %>% unlist() %>% as.numeric()
  
  ## merge/join variants with TSS
  # TSS file: no header, chr, start, stop, gene name, score, strand
  TSS = read.table(file=TSS.file, header=FALSE, sep='\t') %>%
    setNames(c("chr", "start", "stop", "gene", "score", "strand"))
  TSS$center = with(TSS, (start + stop)/2)
  
  # variant file: 1-3 (loc), 4 (hgID), 5 (tissue), 6 (gene symbol), 7 (PIP)
  variants = read.table(file=var.file, header=FALSE, sep='\t') %>%
    setNames(c("chr", "start", "stop", "hgID", "tissue", "gene", "PIP", "TPM"))
  
  # add "center" TSS column and merge to average TSS per gene
  TSS.col = dplyr::select(TSS, c(gene, center)) %>%
    group_by(gene) %>% 
    summarise(center = mean(center))

  var.merged = inner_join(variants, TSS.col, by="gene")

  ## compute TSS-gene distance
  # subtract and absolute value from variant start loc
  var.merged$distance = with(var.merged, abs(center-start))
  
  # add distance group
  var.merged$distanceGroup = 0
  distances = c(0, distances)
  
  # for (i in 1:length(var.merged)){
  #   distance.this = var.merged$distance[i]
  #   for (j in 1:length(distances-2)){
  #     distance.min = distances[i]
  #     distance.max = distances[i+1]
  #     if (distance.this>distance.min & distance.this<=distance.max) {
  #       var.merged$distanceGroup = distance.max
  #     }
  #   }
  # }
  
  for (i in 1:(length(distances)-2)){
    distance.min = distances[i]
    distance.max = distances[i+1]
    var.merged$distanceGroup[var.merged$distance>distance.min & var.merged$distance<=distance.max] = distance.max
  }

  var.merged = dplyr::select(var.merged, -center, -distance)
  
  write.table(var.merged, file=out.file, quote=FALSE, sep='\t', col.names=FALSE, row.names=FALSE)
  

  
}


main()
