suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)})

main <- function() {
  ## get files from snakemake
  countFile = (snakemake@input$countMatrix)
  biosampleFile = (snakemake@input$samples)
  commonVarFile = (snakemake@input$commonVarPerBiosample)
  varPerGTExTissueFile = (snakemake@input$variantsPerGTExTissue)
  GTExTissues = (snakemake@params$GTExTissues) %>% strsplit(" ") %>% unlist %>% sort()
  sampleKeyFile = (snakemake@params$sampleKey)
  outFile = (snakemake@output$enrichmentTable)
  totalCommonVar = 9765167
  
	## read in files
  # count matrix
  biosamples = read.table(biosampleFile, header=FALSE, stringsAsFactors=FALSE)
  # each row = biosample, each col = tissue
	countMatrix = read.table(countFile, header=TRUE, stringsAsFactors=FALSE)

	# common variants by tissue/biosample
	commonVarPerBiosample = read.table(commonVarFile, header=FALSE,stringsAsFactors=FALSE)
	colnames(commonVarPerBiosample) = c('nCommonVariantsOverlappingEnhancers','Biosample')

	# variants by tissue
	variantsByGTExTissue = read.table(varPerGTExTissueFile, header=FALSE, stringsAsFactors=FALSE); 
	colnames(variantsByGTExTissue) = c('GTExTissue','nVariantsGTExTissue')

	# make matrix: columns GTExTissue, Biosample, nVariantsOverlappingEnhancers, nVariantsGTExTissue, nCommonVariantsOverlappingEnhancers,  nCommonVariants, enrichment
	enrMatrix = pivot_longer(countMatrix, cols=-Biosample, names_to='GTExTissue', values_to='nVariantsOverlappingEnhancers')
	enrMatrix[enrMatrix=="Cells_EBV.transformed_lymphocytes"] = "Cells_EBV-transformed_lymphocytes"
	enrMatrix = left_join(enrMatrix, variantsByGTExTissue, by='GTExTissue')
	enrMatrix = left_join(enrMatrix, commonVarPerBiosample, by='Biosample')
	enrMatrix$nCommonVariants = totalCommonVar
	enrMatrix$enrichment = enrMatrix$nVariantsOverlappingEnhancers/enrMatrix$nVariantsGTExTissue/(enrMatrix$nCommonVariantsOverlappingEnhancers/enrMatrix$nCommonVariants)
	
	if (is.na(sampleKeyFile) || sampleKeyFile=="None"){
	  enrMatrix$sampleName = enrMatrix$Biosample
	} else {
	  cat.data = read.table(sampleKeyFile, header=TRUE, sep="\t", fill=TRUE)
	  if ("sampleName" %in% colnames(cat.data)){
	    cat.data = dplyr::select(cat.data, biosample, sampleName)
	    enrMatrix = left_join(enrMatrix, cat.data, by=c('Biosample'='biosample'))
	  } else {
	    enrMatrix$sampleName = enrMatrix$Biosample
	  }
	}
	
	# write table
	write.table(enrMatrix, file=outFile, sep="\t", quote=F, row.names=F, col.names=T)
}

main()
