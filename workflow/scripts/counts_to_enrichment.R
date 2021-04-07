# input arguments: count matrix, common var per  biosample, variants per GTEx tissue, total common variants, variantSet
# output matrix with columns: GTExTissue, Biosample, nVariantsOverlappingEnhancers, nVariantsGTExTissue, nCommonVariantsOverlappingEnhancers, nCommonVariants, enrichment

library(dplyr); library(tidyr); library(optparse)
  
  main <- function() {
    option_list <- list(
      make_option(c("--counts"), type="character", default=NA, help="count matrix of variants in each GTEx tissue/biosample overlap"),
      make_option(c("--commonvar"), type="character", default=NA, help="number of common variants overlapping enhancers in each biosample"),
      make_option(c("--varbytissue"), type="character", default=NA, help="number of GTEx variants per tissue"),
      make_option(c("--totalcommonvar"), type="numeric", default=NA, help="total common variants"))
    opt = parse_args(OptionParser(option_list=option_list))
    
	# read in files
	countMatrix = read.csv(opt$counts, sep = '\t', header=FALSE,stringsAsFactors = FALSE) %>% dplyr::select(-"V34")
	commonVarByBiosample = read.csv(opt$commonvar, sep = '\t', header=FALSE,stringsAsFactors = FALSE); 
	colnames(commonVarByBiosample) = c('Biosample','nCommonVariantsOverlappingEnhancers')
	variantsByGTExTissue = read.csv(opt$varbytissue, sep = '\t', header=FALSE,stringsAsFactors = FALSE); 
	colnames(variantsByGTExTissue) = c('GTExTissue','nVariantsGTExTissue')
	totalCommonVar = opt$totalcommonvar

	# make matrix: columns GTExTissue, Biosample, nVariantsOverlappingEnhancers, nVariantsGTExTissue, nCommonVariantsOverlappingEnhancers,  nCommonVariants, enrichment
	colnames(countMatrix) = c('Biosample', sort(variantsByGTExTissue$GTExTissue))

	enrMatrix = gather(countMatrix, key='GTExTissue', value='nVariantsOverlappingEnhancers', -Biosample) %>% left_join(variantsByGTExTissue) %>% left_join(commonVarByBiosample)
	enrMatrix$nCommonVariants = totalCommonVar
	enrMatrix$enrichment = enrMatrix$nVariantsOverlappingEnhancers/enrMatrix$nVariantsGTExTissue/(enrMatrix$nCommonVariantsOverlappingEnhancers/enrMatrix$nCommonVariants)

	# write table
	write.table(enrMatrix, file="", sep="\t", quote=F, row.names=F, col.names=T)
}

main()
