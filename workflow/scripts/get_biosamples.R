suppressPackageStartupMessages({library(dplyr)
  library(tidyr)
  library(stringr)})

  main <- function() {
  ## get files from snakemake
  biosampleKey = (snakemake@input$biosampleKey)
  outFile = (snakemake@output$samples)

	key = read.table(file=biosampleKey, sep="\t", fill=TRUE, header=TRUE)
	biosamples = key$biosample[order(key$biosample)] # alphabetized?
	print(biosamples)
   
	# write matrix
	write.table(biosamples, file=outFile, sep="\t", quote=F, row.names=F, col.names=F)
}

main()
