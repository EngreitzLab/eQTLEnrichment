suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)})


main <- function() {
  var.file = (snakemake@input$filteredGTExVariants) 
  distances = (snakemake@params$distances) %>% unlist() %>% as.numeric()
  GTExTissue = (snakemake@params$GTExTissues) %>% strsplit(" ") %>% unlist()
  out.file.all = (snakemake@output$variantsPerGTExTissue)
  out.files.by.dist = (snakemake@params$variantsPerGTExTissueByDist) %>% strsplit(" ") %>% unlist()
  
  variants = read.table(file=var.file, header=FALSE, sep='\t') %>%
    setNames(c("chr", "start", "stop", "hgID", "tissue", "gene", "PIP", "TPM", "distance")) 
  
  # df for all variants
  df.all = data.frame(GTExTissue)
  df.all$nVariants = 0
  for (i in length(GTExTissue)){
    GTExTissue.this = GTExTissue[i]
    variants.this = dplyr::filter(variants, tissue==GTExTissue.this) %>% dplyr::select(chr, start, stop) %>% distinct() %>% nrow()
    df.all$nVariants[i] = variants.this
  }
  write.table(df.all, file=out.file.all, sep="\t", row.names=TRUE, quote=FALSE)
  
  # by distance
  for (i in length(distances)){
    distance.this = distances[i]
    variants.these = dplyr::filter(variants, distance<=distance.this)
    df.this = data.frame(GTExTissue)
    df.this$nVariants = 0
    for (k in length(GTExTissue)){
      GTEx.Tissue.this = GTExTissue[k]
      variants.this = dplyr::filter(variants.these, tissue==GTExTissue.this) %>% dplyr::select(chr, start, stop) %>% distinct() %>% nrow()
      df.this$nVariants[k] = variants.this
    }
    write.table(df.this, file=out.files.by.dist[i], sep="\t", row.names=TRUE, quote=FALSE)
  }
  
  
  
  
  
# df per distance
  
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
  out.this = dplyr::select(filtered.this, -c(distance, center)) %>% unique()
  write.table(out.this, file=out.files[1], quote=FALSE, sep='\t', col.names=FALSE, row.names=FALSE)
  
  for (i in 2:length(distances)) {
    if (distances[i]==30e6){ # this means use all variants
      filtered.this = var.merged
      out.this = dplyr::select(filtered.this, -c(distance, center)) %>% unique()
      write.table(out.this, file=out.files[i], quote=FALSE, sep='\t', col.names=FALSE, row.names=FALSE)
    } else {
      filtered.this = dplyr::filter(var.merged, var.merged$distance<=distances[i], var.merged$distance>=distances[i-1])
      out.this = dplyr::select(filtered.this, -c(distance, center)) %>% unique()
      write.table(out.this, file=out.files[i], quote=FALSE, sep='\t', col.names=FALSE, row.names=FALSE)
    }
    
  }

  
}


main()
