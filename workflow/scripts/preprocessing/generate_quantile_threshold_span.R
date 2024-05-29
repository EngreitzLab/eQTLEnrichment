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
nSteps = snakemake@params$nSteps %>% as.numeric()
binary = snakemake@params$binary
biosamples = snakemake@params$biosamples %>% as.character() %>% strsplit(" ") %>% unlist() 
varInt_files = snakemake@input$varInt %>% as.character() %>% strsplit(" ") %>% unlist() 

outFile = snakemake@output$outFile

if (binary %in% c("TRUE", "True", TRUE)){
  print("binary")
  nSteps=2
} else {
  print("non-binary")
}

sampleKey$GTExTissue = sampleKey$tissue
sampleKey$key = paste0(sampleKey$GTExTissue, ".", sampleKey$biosample)

col_names = c("varChr", "varStart", "varEnd", "variantID", "eGene", "GTExTissue", "PIP", "distanceBin", "enhChr", "enhStart", "enhEnd", "biosample", "TargetGene", "score" )

for (i in 1:nrow(sampleKey)){
	GTExTissue_this = sampleKey$GTExTissue[i]
	biosample_this = sampleKey$biosample[i]
	biosample_idx = which(biosamples==biosample_this)
	varInt_this = varInt_files[biosample_idx]
	varInt = read.table(file=varInt_this, sep="\t", header=FALSE) %>% setNames(col_names)
	varInt_filtered = dplyr::filter(varInt, GTExTissue==GTExTissue_this, eGene==TargetGene) %>%
  		dplyr::select(varChr, varStart, varEnd, eGene, score) %>% distinct()
	if (i==1){varInt_all = varInt_filtered} else {
		varInt_all = rbind(varInt_all, varInt_filtered)
	}
}

scores = as.numeric(varInt_all$score)
prob_vector = seq(0, 1, length.out=nSteps)
thresholds = data.frame(unname(quantile(scores, na.rm=T, probs=prob_vector))) %>% distinct()

colnames(thresholds) = "threshold"

write.table(thresholds, outFile, row.names=FALSE, col.names=FALSE, sep="\t")
