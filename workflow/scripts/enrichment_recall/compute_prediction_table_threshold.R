suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(data.table)})

varPredInt_files = (snakemake@input$varPredInt) %>% strsplit(" ") %>% unlist() 
GTExVariantsFile = (snakemake@input$filteredGTExVariantsFinal)
thresholdSpan = read.table(snakemake@input$thresholdSpan, sep="\t", header=FALSE) %>% setNames("threshold")
biosamples = (snakemake@params$biosamples) %>% strsplit(" ") %>% unlist()
map = fread(snakemake@input$map, sep="\t", header=TRUE) #tissue, biosample
outFile = (snakemake@output$predTable)
method.this = snakemake@wildcards$method

# get number of variant-gene pairs per tissue 
GTExVariants = fread(file=GTExVariantsFile, sep="\t", header=FALSE) %>%
  setNames(c("varChr", "varStart", "varEnd", "variantID",  "eGene", "GTExTissue", "PIP", "distanceBin"))
tissues = unique(GTExVariants$GTExTissue)
nVariantsTotal = data.frame(GTExTissue=tissues)
nVariantsTotal$total.variants = 0
for (i in 1:length(tissues)){
	nVariantsTotal$total.variants[i] = dplyr::filter(GTExVariants, GTExTissue==tissues[i]) %>% dplyr::select(varChr, varStart, varEnd, eGene) %>% distinct() %>% nrow()
}

# get number of correctly-predicted variants
for (b in 1:length(biosamples)){ # iterate through biosamples
	sample.this = biosamples[b]
    varIntFile = varPredInt_files[b]
    
    size.file = file.info(varIntFile)$size
    size.threshold = 100 # in bytes, for empty file
	
	counts_overlap = data.frame(matrix(ncol=length(tissues), nrow=length(thresholdSpan$threshold))) # number of eVariant-eGene pairs where eVariant overlaps predicted enhancer
	counts_linking = data.frame(matrix(ncol=length(tissues), nrow=length(thresholdSpan$threshold))) # number of eVariant-eGene pairs where eVariant overlaps predicted enhancer linked to correct eGene
	colnames(counts_overlap) = tissues; colnames(counts_linking) = tissues
	counts_overlap$Biosample = sample.this; counts_linking$Biosample = sample.this
	counts_overlap$threshold = thresholdSpan$threshold; counts_linking$threshold = thresholdSpan$threshold

	if (size.file<size.threshold){ # file is empty (no intersection)
     		 for (t in tissues)
			 	counts_overlap[[t]] = 0; counts_linking[[t]] = 0
    } else { # read in file
		variantsInt = fread(varIntFile, header=FALSE, sep="\t") %>% 
						setNames(c("varChr", "varStart", "varEnd", "variantID", "eGene", "GTExTissue", "PIP", "distance_bin",
							"enhChr", "enhStart", "enhEnd", "Biosample", "TargetGene", "score"))
		for (n in 1:nrow(thresholdSpan)){ 	# iterate through thresholds
			score.thresh = thresholdSpan$threshold[n]
			variantsThresh = dplyr::filter(variantsInt, score>=score.thresh)  %>% dplyr::select(varChr, varStart, varEnd, GTExTissue, eGene, TargetGene)
			for (t in tissues){
				variantsTissue = dplyr::filter(variantsThresh, GTExTissue==t)
				overlap =  dplyr::select(variantsTissue, varChr, varStart, varEnd, eGene) %>% distinct() %>% nrow()
				linking = dplyr::filter(variantsTissue, eGene==TargetGene) %>% dplyr::select(varChr, varStart, varEnd, eGene) %>% distinct() %>% nrow()

				counts_overlap[counts_overlap$threshold==score.thresh,t] = overlap
				counts_linking[counts_overlap$threshold==score.thresh,t] = linking
      		}
    	}
		rm(variantsInt)
	}
	if (b==1){
		counts_overlap_all = counts_overlap; counts_linking_all = counts_linking
	} else {
		counts_overlap_all = rbind(counts_overlap_all, counts_overlap); counts_linking_all = rbind(counts_linking_all, counts_linking)
  }
}

df_overlap = pivot_longer(counts_overlap_all, cols=-c(Biosample,threshold), names_to='GTExTissue', values_to='nVariantGenePairsOverlappingEnhancers')
df_linking = pivot_longer(counts_linking_all, cols=-c(Biosample,threshold), names_to='GTExTissue', values_to='nVariantsOverlappingEnhancersCorrectGene')
df = left_join(df_overlap, df_linking, by=c("Biosample", "GTExTissue", "threshold")) %>% 
	left_join(nVariantsTotal, by="GTExTissue")

# calculate metrics across all matches
map$key = paste0(map$tissue, ".", map$biosample)

df_all = dplyr::mutate(df, key=paste0(GTExTissue, ".", Biosample)) %>%
	dplyr::filter(key %in% map$key) %>%
	group_by(threshold) %>%
	summarize(nVariantGenePairsOverlappingEnhancers = sum(nVariantGenePairsOverlappingEnhancers), nVariantsOverlappingEnhancersCorrectGene = sum(nVariantsOverlappingEnhancersCorrectGene), total.variants = sum(total.variants))  %>%
	mutate(GTExTissue="all_matches", Biosample="all_matches")
df = rbind(df, df_all)
	
# calculate
df$recall.total = with(df, nVariantGenePairsOverlappingEnhancers/total.variants) # fraction of variants overlapping enhancers
df$recall.linking = with(df, nVariantsOverlappingEnhancersCorrectGene/total.variants)  # fraction of variants overlapping enhancers linked to correct gene
df$correctGene.ifOverlap = with(df, nVariantsOverlappingEnhancersCorrectGene/nVariantGenePairsOverlappingEnhancers) # fraction of variants linked to correct gene GIVEN they overlap enhancers 
df$method = method.this

fwrite(df, file=outFile, sep="\t", quote=F, row.names=F, col.names=T)

