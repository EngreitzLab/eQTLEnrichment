suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(data.table)})

## get files from snakemake
method = (snakemake@wildcards$method)
thresholdSpan = (read.table(snakemake@input$thresholdSpan, sep='\t', header=FALSE)) %>% setNames("threshold")
countFile = (snakemake@input$countMatrix)
biosamples = (snakemake@params$biosamples) %>% strsplit(" ") %>% unlist()
sign_threshold = snakemake@params$thresholdPval %>% as.numeric()
map = (snakemake@input$map)
varPerGTExTissueFile = (snakemake@input$variantsPerGTExTissue)
commonVarInt_files = (snakemake@input$commonVarInt)
outFile = (snakemake@output$enrichmentTable)

bgCommonVar_file = (snakemake@input$commonVarCount)
bgCommonVar = read.table(bgCommonVar_file, header=FALSE, sep="\t") %>% setNames("n")
bgCommonVar = as.numeric(bgCommonVar$n)

# count matrix: each row = biosample, each col = tissue
countMatrix = read.table(file=countFile, header=TRUE, stringsAsFactors=FALSE)
countMatrix_pivot = pivot_longer(countMatrix, cols=-c(Biosample,threshold), names_to='GTExTissue', values_to='nVariantsOverlappingEnhancers')

# variants per GTEx tissue
variantsByGTExTissue = read.table(varPerGTExTissueFile, header=TRUE, stringsAsFactors=FALSE, sep="\t") # col names: tissue, n_all, n_bin1, etc.
variantsByGTExTissue = dplyr::select(variantsByGTExTissue, tissue, n_all)
colnames(variantsByGTExTissue) = c("GTExTissue", "nVariantsGTExTissue")

for (b in 1:length(biosamples)){
	sample.this = biosamples[b]
	commonVarPredIntFile = commonVarInt_files[b]
	commonVarPredInt = read.table(file=commonVarPredIntFile, header=TRUE, stringsAsFactors=FALSE) %>%
      setNames(c("varChr", "varStart", "varEnd", "rsID", "enhChr", "enhStart", "enhEnd", "Biosample", "TargetGene", "score"))

	commonVarPerBiosample =data.frame(threshold = thresholdSpan$threshold)
	commonVarPerBiosample$Biosample = sample.this
	commonVarPerBiosample$nCommonVariantsOverlappingEnhancers = 0

	for (t in 1:nrow(thresholdSpan)){
		score.thresh = thresholdSpan$threshold[t]
		commonVar_this = dplyr::filter(commonVarPredInt, score>= score.thresh)
		count = dplyr::select(commonVar_this, rsID) %>% distinct() %>% nrow()
		commonVarPerBiosample$nCommonVariantsOverlappingEnhancers[t] = count
	}

	# enrichment matrix for this biosample
	enrMatrix = dplyr::filter(countMatrix_pivot, Biosample==sample.this)
	enrMatrix[enrMatrix=="Cells_EBV.transformed_lymphocytes"] = "Cells_EBV-transformed_lymphocytes"
	enrMatrix = left_join(enrMatrix, variantsByGTExTissue, by='GTExTissue')
	enrMatrix = left_join(enrMatrix, commonVarPerBiosample, by=c("Biosample", "threshold"))
	enrMatrix$nCommonVariants = bgCommonVar
	enrMatrix$enrichment = enrMatrix$nVariantsOverlappingEnhancers/enrMatrix$nVariantsGTExTissue/(enrMatrix$nCommonVariantsOverlappingEnhancers/enrMatrix$nCommonVariants)

	if (b==1){
    enrMatrix_all = enrMatrix
	} else {
    enrMatrix_all = rbind(enrMatrix_all, enrMatrix)
  	}
}

## compute aggregate enrichment for matching tissues/biosamples across thresholds
map = fread(map, sep="\t", header=TRUE)
map$key = paste0(map$tissue, ".", map$biosample)

summ = enrMatrix_all
summ$match_id = paste0(summ$GTExTissue, ".", summ$Biosample)
summ = dplyr::filter(summ, match_id %in% map$key) %>%
	dplyr::select(threshold, nVariantsGTExTissue, nVariantsOverlappingEnhancers, nCommonVariants, nCommonVariantsOverlappingEnhancers) %>%
	group_by(threshold) %>%
	summarize(nVarGTExTissue_sum = sum(nVariantsGTExTissue), 
						nVariantsOverlappingEnhancers_sum = sum(nVariantsOverlappingEnhancers),
						nCommonVariants_sum = sum(nCommonVariants),
						nCommonVariantsOverlappingEnhancers_sum = sum(nCommonVariantsOverlappingEnhancers)) %>%
	mutate(enrichment = nVariantsOverlappingEnhancers_sum/nVarGTExTissue_sum/(nCommonVariantsOverlappingEnhancers_sum/nCommonVariants_sum),
				GTExTissue = "all_matches",
				Biosample = "all_matches",
				nVariantsGTExTissue = nVarGTExTissue_sum,
				nVariantsOverlappingEnhancers = nVariantsOverlappingEnhancers_sum,
				nCommonVariantsOverlappingEnhancers = nCommonVariantsOverlappingEnhancers_sum,
				nVariantsGTExTissue = nVarGTExTissue_sum,
				nCommonVariants = nCommonVariants_sum) %>%
	dplyr::select(-c(nVarGTExTissue_sum, nVariantsOverlappingEnhancers_sum, nCommonVariantsOverlappingEnhancers_sum, nVarGTExTissue_sum,nCommonVariants_sum))

enrMatrix_all = rbind(enrMatrix_all, summ)


## stats about risk ratio (RR) aka enrichment
# calculate CI of RR and SE(log RR); see: https://sphweb.bumc.bu.edu/otlt/mph-modules/bs/bs704_confidence_intervals/bs704_confidence_intervals8.html
z = qnorm(sign_threshold/2, lower.tail=FALSE) # e.g. 1.96 for p=0.05
calcs = enrMatrix_all
calcs$n1 = calcs$nVariantsGTExTissue
calcs$x1 = calcs$nVariantsOverlappingEnhancers
calcs$n2 = calcs$nCommonVariants
calcs$x2 = calcs$nCommonVariantsOverlappingEnhancers

calcs$SE_log_enr = with(calcs, sqrt(((n1-x1)/x1)/n1) + ((n2-x2)/x2)/n2)
calcs$log_enr = log(calcs$enrichment)
calcs$log_CI_enr_low = with(calcs, log_enr - z*SE_log_enr)
calcs$log_CI_enr_high = with(calcs, log_enr + z*SE_log_enr)
calcs$CI_enr_low = exp(calcs$log_CI_enr_low)
calcs$CI_enr_high = exp(calcs$log_CI_enr_high)

# significance (from Jesse's CredibleSetTools.R)
calcs$p = with(calcs, mapply(FUN=phyper, x1, n1, n2, x1+x2, log.p=FALSE, lower.tail=FALSE))
calcs$p.adjust = p.adjust(calcs$p, method="bonferroni")

enrMatrix_all$CI_enr_low = calcs$CI_enr_low
enrMatrix_all$CI_enr_high = calcs$CI_enr_high
enrMatrix_all$SE_log_enr = calcs$SE_log_enr
enrMatrix_all$p_adjust_enr = calcs$p.adjust
  
## write output
fwrite(enrMatrix_all, file=outFile, sep="\t", quote=F, row.names=F, col.names=T)
