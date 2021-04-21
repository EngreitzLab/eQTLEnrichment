suppressPackageStartupMessages({library(dplyr)
  library(tidyr)
  library(stringr)})

  main <- function() {
  ## get files from snakemake
  dfFile = (snakemake@input$variantsPredictionsInt)
  biosampleFile = (snakemake@input$samples)
  GTExTissues = (snakemake@params$GTExTissues) %>% strsplit(" ") %>% unlist %>% sort()
  outFile = (snakemake@output$countMatrix)

	## read in files
  biosamples = read.table(biosampleFile, header=FALSE, stringsAsFactors=FALSE); biosamples = biosamples[[1]]
  df = read.table(gzfile(dfFile), header=FALSE)
  colnames(df) = c('var.chr','var.start','var.end','hgid','GTExTissue','eGene','PIP', 'TPM', 
                   'enh.chr', 'enh.start', 'enh.end', 'Biosample', 'TargetGene')
  df = dplyr::select(df, hgid, GTExTissue, Biosample)
  
  # initialize counts matrix
  counts = data.frame(matrix(ncol=length(GTExTissues), nrow=length(biosamples)))
  colnames(counts)=GTExTissues
  counts$Biosample=biosamples
  
  # get counts
  for (tissue in GTExTissues){
    vars.tissue=filter(df, GTExTissue==tissue)
    for (sample in biosamples){
      vars.int = filter(vars.tissue, Biosample==sample) %>% distinct()
      counts[counts$Biosample==sample,tissue] = nrow(vars.int)
    }
  }
  
	# write matrix
	write.table(counts, file=outFile, sep="\t", quote=F, row.names=F, col.names=T)
}

main()
