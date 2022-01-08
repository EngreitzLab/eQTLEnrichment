suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)})

main <- function() {
  ## get files from snakemake
  method = (snakemake@wildcards$method)
  threshold = (snakemake@wildcards$threshold)
  outDir = (snakemake@params$outDir)
  countFile = (snakemake@input$countMatrix)
  biosampleFile = (snakemake@input$samples)
  varPerGTExTissueFile = (snakemake@input$variantsPerGTExTissue)
  GTExTissues = (snakemake@params$GTExTissues) %>% strsplit(" ") %>% unlist %>% sort()
  sampleKeyFile = (snakemake@params$sampleKey)
  outFile = (snakemake@output$enrichmentTable)
  totalCommonVar = 9765167
  
	## read in files
  # count matrix
  biosamples = read.table(file=biosampleFile, header=FALSE, stringsAsFactors=FALSE) %>% setNames(c("Biosample"))

  # each row = biosample, each col = tissue
	countMatrix = read.table(file=countFile, header=TRUE, stringsAsFactors=FALSE)

	# common variants by tissue/biosample
	# loop through biosamples and make this 
	for (i in 1:nrow(biosamples)){
    sample.this = biosamples$Biosample[i]
    # read in count file
    file.this = file.path(outDir, method, sample.this, 
                          paste0("commonVarPerBiosample.", threshold, ".tsv"))
    file.size = file.info(file.this)$size
    size.threshold = 10
    if (file.size<size.threshold){
      counts.this = data.frame(matrix(nrow=1, ncol=2)) %>% setNames(c('nCommonVariantsOverlappingEnhancers','Biosample'))
      counts.this$nCommonVariantsOverlappingEnhancers[1]=0
      counts.this$Biosample[1]=sample.this
    } else {
      counts.this = read_table(file=file.this, col_names=FALSE) %>% setNames(c('nCommonVariantsOverlappingEnhancers','Biosample'))
    }
    
     if (i==1){
      commonVarPerBiosample=counts.this
    }
    else{
      commonVarPerBiosample=rbind(commonVarPerBiosample, counts.this)
    }
	}

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
