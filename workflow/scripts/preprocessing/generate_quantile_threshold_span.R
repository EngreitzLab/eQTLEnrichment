suppressPackageStartupMessages({
  library(plyr)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(forcats)
  library(data.table)
  library(stringr)
})

### INPUTS
method_this = snakemake@wildcards$method
sampleKey = fread(file=snakemake@input$map, sep="\t")
n_steps_target = snakemake@params$nSteps %>% as.numeric()
binary = snakemake@params$binary
provided_threshold = snakemake@params$threshold %>% as.numeric()
biosamples = snakemake@params$biosamples %>% as.character() %>% strsplit(" ") %>% unlist() 
varInt_files = snakemake@input$varInt %>% as.character() %>% strsplit(" ") %>% unlist() 

outFile = snakemake@output$outFile

if (binary %in% c("TRUE", "True", TRUE)){
	print("binary")
	nSteps=2
	thresholds_all = data.frame(threshold = c(0, 1))
} else {
	print("non-binary")
	sampleKey$GTExTissue = sampleKey$tissue
	sampleKey$key = paste0(sampleKey$GTExTissue, ".", sampleKey$biosample)

	col_names = c("varChr", "varStart", "varEnd", "variantID", "eGene", "GTExTissue", "PIP", "distanceBin", "enhChr", "enhStart", "enhEnd", "biosample", "TargetGene", "score" )

	varInt_list = vector(mode = "list", length = nrow(sampleKey))
	for (i in 1:nrow(sampleKey)){
		GTExTissue_this = sampleKey$GTExTissue[i]
		biosample_this = sampleKey$biosample[i]
		biosample_idx = which(biosamples==biosample_this)
		varInt_this = varInt_files[biosample_idx]
		if (file.size(varInt_this) > 100) {
			varInt = fread(file=varInt_this, sep="\t", header=FALSE) %>% setNames(col_names)
			varInt_filtered = dplyr::filter(varInt, GTExTissue==GTExTissue_this, eGene==TargetGene) %>%
				dplyr::select(varChr, varStart, varEnd, eGene, score) %>% 
				dplyr::filter(!is.na(score), is.finite(score)) %>%
				distinct()
			varInt_list[[i]] = varInt_filtered
		}	
	}
	print(varInt_list)
	varInt_all = rbindlist(varInt_list) %>% as_tibble()

	# evenly-spaced scores (n_steps_target)
	scores = as.numeric(varInt_all$score)
	n_distinct_scores = unique(scores) %>% length()
	even_spaced_scores = data.frame(threshold = seq(min(scores), max(scores), length.out=n_steps_target))

	# quantile spaced scores (n_steps_target)
	prob_vector = seq(0, 1, length.out=n_steps_target)
	quantile_scores = unname(quantile(scores, na.rm=T, probs=prob_vector)) 
	quantile_scores = c(quantile_scores, provided_threshold)

	thresholds = data.frame(threshold=quantile_scores) %>% distinct()

	n_steps = n_steps_target
	while ((nrow(thresholds)<n_steps_target) & (n_steps<n_distinct_scores)){ # less than target # thresholds, but not max-ed out
		n_steps = n_steps + 1
		prob_vector = seq(0, 1, length.out = n_steps)
		quantile_scores = unname(quantile(scores, na.rm=T, probs=prob_vector))
		quantile_scores = c(quantile_scores, provided_threshold)

		thresholds = data.frame(threshold=quantile_scores) %>% distinct()
	}

	# combine
	thresholds_all = rbind(thresholds, even_spaced_scores) %>%
		distinct() %>%
		arrange(threshold)
}


print(thresholds_all)

write.table(thresholds_all, outFile, row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
