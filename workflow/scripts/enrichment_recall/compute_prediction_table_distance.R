suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(data.table)})

varPredInt_files = (snakemake@input$varPredInt) %>% strsplit(" ") %>% unlist() 
GTExVariantsFile = (snakemake@input$filteredGTExVariantsFinal)
score.thresh = (snakemake@params$threshold) %>% as.numeric()
distances_min = (snakemake@params$distances_min) %>% as.character() %>% strsplit(" ") %>% unlist() %>% as.numeric()
distances_max = (snakemake@params$distances_max) %>% as.character() %>% strsplit(" ") %>% unlist() %>% as.numeric()
outFile = (snakemake@output$predTable)
biosamples = (snakemake@params$biosamples) %>% strsplit(" ") %>% unlist()
method.this = snakemake@wildcards$method

# get number of variants per tissue
GTExVariants = fread(file=GTExVariantsFile, sep="\t", header=FALSE) %>%
  setNames(c("varChr", "varStart", "varEnd", "variantID",  "eGene", "GTExTissue", "PIP", "distanceBin"))
tissues = unique(GTExVariants$GTExTissue)
for (d in 1:length(distances_max)){
	if (distances_max[d]==30000) {
		variants.dist = GTExVariants
	} else {
		variants.dist = dplyr::filter(GTExVariants, distanceBin == distances_max[d])
	}
	nVar.this = dplyr::select(variants.dist, varChr, varStart, varEnd, eGene, GTExTissue) %>% distinct() %>%
		group_by(GTExTissue) %>% summarize(total.variants = n())
	nVar.this$distance_min = distances_min[d]
	nVar.this$distance_max = distances_max[d]

	if (d==1) {nVariantsTotal  = nVar.this} else {
		nVariantsTotal = rbind(nVariantsTotal, nVar.this)
	}
}

for (d in 1:length(distances_max)){
	for (b in 1:length(biosamples)){ # iterate through biosamples
		sample.this = biosamples[b]
		varIntFile = varPredInt_files[b]
		
		size.file = file.info(varIntFile)$size
		size.threshold = 100 # in bytes, for empty file
		
		df = data.frame(GTExTissue = tissues)
		df$nVariantGenePairsOverlappingEnhancers = 0
		df$nVariantsOverlappingEnhancersCorrectGene = 0
		df$Biosample = sample.this
		df$distance_min = distances_min[d]
		df$distance_max = distances_max[d]

		if (size.file>size.threshold){ # file is not empty -> read in
			variantsInt = fread(varIntFile, header=FALSE, sep="\t") %>% 
							setNames(c("varChr", "varStart", "varEnd", "variantID", "eGene", "GTExTissue", "PIP", "distance_bin",
								"enhChr", "enhStart", "enhEnd", "Biosample", "TargetGene", "score"))
			variantsInt.thresh = dplyr::filter(variantsInt, score>=score.thresh)
			rm(variantsInt)
			if (distances_max[d]==30000){
				variantsInt.dist = variantsInt.thresh
			} else {
				variantsInt.dist = dplyr::filter(variantsInt.thresh, distance_bin==distances_max[d])
			}

			for (t in 1:length(tissues)){ # iterate through tissues
				variantsInt.tissue = dplyr::filter(variantsInt.dist, GTExTissue==tissues[t])
				df$nVariantGenePairsOverlappingEnhancers[t] =  dplyr::select(variantsInt.tissue, varChr, varStart, varEnd, eGene) %>% distinct() %>% nrow()
				df$nVariantsOverlappingEnhancersCorrectGene[t] = dplyr::filter(variantsInt.tissue, eGene==TargetGene) %>% dplyr::select(varChr, varStart, varEnd, eGene) %>% distinct() %>% nrow()
			}
		}
		# concat across bisamples
		if (b==1) {df_biosample = df} else {df_biosample = rbind(df_biosample, df)}
	}
	#concat across distances
	if (d==1) {df_all = df_biosample} else {df_all = rbind(df_all, df_biosample)}
}

df_all = left_join(df_all, nVariantsTotal, by=c("distance_min", "distance_max", "GTExTissue")) # add total.variants

# calculate
df_all$recall.total = with(df_all, nVariantGenePairsOverlappingEnhancers/total.variants) # fraction of variants overlapping enhancers
df_all$recall.linking = with(df_all, nVariantsOverlappingEnhancersCorrectGene/total.variants)  # fraction of variants overlapping enhancers linked to correct gene
df_all$correctGene.ifOverlap = with(df_all, nVariantsOverlappingEnhancersCorrectGene/nVariantGenePairsOverlappingEnhancers) # fraction of variants linked to correct gene GIVEN they overlap enhancers 

df_all$method = method.this

fwrite(df_all, file=outFile, sep="\t", quote=F, row.names=F, col.names=T)

