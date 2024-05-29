suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(data.table)})


main <- function() {
  var.file = (snakemake@input$filteredGTExVariants) 
  TSS.file = (snakemake@params$TSS)
  out.file = (snakemake@output$GTExVariantsDistance)
  distances = (snakemake@params$distances) %>% as.character() %>% strsplit(" ") %>% unlist() %>% as.numeric()
  
  ## merge/join variants with TSS
  # TSS file: no header, chr, start, stop, gene name, score, strand
  TSS = read.table(file=TSS.file, header=FALSE, sep='\t') %>%
    setNames(c("chr", "start", "end", "gene", "score", "strand"))
  TSS$center = with(TSS, (start + end)/2)

  # variant file:chr,start,end,variant_id,gene_hgnc,tissue,pip
  variants = read.table(file=var.file, header=FALSE, sep='\t') %>%
    setNames(c("chr", "start", "end", "variant_id", "gene", "tissue", "pip"))

  # add "center" TSS column and merge to average TSS per gene
  TSS.col = dplyr::select(TSS, c(gene, center)) %>%
    group_by(gene) %>% 
    summarise(center = mean(center))

  var.merged = inner_join(variants, TSS.col, by="gene")
print(head(var.merged))
  ## compute TSS-gene distance
  # subtract and absolute value from variant start loc
  var.merged$distance = with(var.merged, abs(center-start))
  
  # add distance group
  var.merged$distanceGroup = 0
  distances = c(0, distances, 30000) # want bins: 0-a, a-b, b-c, over c (n+1 categories)
  
  
  for (i in 1:(length(distances)-1)){
    distance.min = distances[i] * 1000 
    distance.max = distances[i+1] * 1000
    var.merged$distanceGroup[var.merged$distance>distance.min & var.merged$distance<=distance.max] = distance.max/1000
  }

  var.merged = dplyr::select(var.merged, -center, -distance)
  fwrite(var.merged, file=out.file, quote=FALSE, sep='\t', col.names=FALSE, row.names=FALSE)

  
}


main()
