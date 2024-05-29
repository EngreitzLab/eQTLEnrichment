suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)})

var.file = (snakemake@input$filteredGTExVariantsFinal) 
distances_max = (snakemake@params$distances_max) %>% as.character() %>% strsplit(" ") %>% unlist() %>% as.numeric()
out.file.all = (snakemake@output$variantsPerTissue)

variants = read.table(file=var.file, header=FALSE, sep='\t') %>%
  setNames(c("chr", "start", "end", "variant_id", "gene", "tissue", "PIP", "distance_bin")) 

# calculate n variants in diff categories
all_tissues = unique(variants$tissue)
df = data.frame(tissue=all_tissues)

# initialize count columns for each distance bin
df$n_all = 0

for (d in 1:length(distances_max)) {
	col_name = paste0("n_bin", d)
	df[[col_name]] = 0
}

# iterate through tissues
for (i in 1:nrow(df)) {
  GTExTissue.this = df$tissue[i]
  variants.tissue = dplyr::filter(variants, tissue==GTExTissue.this)
  count = dplyr::select(variants.tissue, chr, start, end) %>% distinct() %>% nrow()
  df$n_all[i] = count

  # iterate through distance bins
  for (d in 1:length(distances_max)){
	col_name = paste0("n_bin", d)
	variants.this = dplyr::filter(variants.tissue, distance_bin==distances_max[d])
	count = dplyr::select(variants.this, chr, start, end) %>% distinct() %>% nrow()

	df[[col_name]][i] = count
  }
}
write.table(df, file=out.file.all, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)