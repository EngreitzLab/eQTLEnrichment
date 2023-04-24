suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)})

var.file = (snakemake@input$filteredGTExVariants) 
distances = (snakemake@params$distances) %>% as.character() %>% strsplit(" ") %>% unlist() %>% as.numeric()
GTExTissue = (snakemake@params$GTExTissues)
out.file.all = (snakemake@output$variantsPerGTExTissue)
out.files.by.dist = (snakemake@output$variantsPerGTExTissueByDist) %>% strsplit(" ") %>% unlist()

variants = read.table(file=var.file, header=FALSE, sep='\t') %>%
  setNames(c("chr", "start", "stop", "hgID", "tissue", "gene", "PIP", "TPM", "distance")) 

# df for all variants
df.all = data.frame(GTExTissue)
df.all$nVariants = 0
for (i in 1:length(GTExTissue)){
  GTExTissue.this = GTExTissue[i]
  variants.this = dplyr::filter(variants, tissue==GTExTissue.this) %>% dplyr::select(chr, start, stop) %>% distinct() %>% nrow()
  df.all$nVariants[i] = variants.this
}
write.table(df.all, file=out.file.all, sep="\t", row.names=TRUE, quote=FALSE)

# by distance

for (i in 1:length(distances)){
  distance.this = distances[i]
  if (distance.this == 30000000) {
    variants.these = variants
  } else {
    variants.these = dplyr::filter(variants, distance==distance.this)
  }
  df.this = data.frame(GTExTissue)
  df.this$nVariants = 0
  for (k in 1:length(GTExTissue)){
    GTExTissue.this = GTExTissue[k]
    variants.this = dplyr::filter(variants.these, tissue==GTExTissue.this) %>% dplyr::select(chr, start, stop) %>% distinct() %>% nrow()
    df.this$nVariants[k] = variants.this
  }

  write.table(df.this, file=out.files.by.dist[i], sep="\t", row.names=TRUE, quote=FALSE)
}
