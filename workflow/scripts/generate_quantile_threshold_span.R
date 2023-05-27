suppressPackageStartupMessages({
  library(plyr)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(forcats)
  library(data.table)
})

### INPUTS
method_this = snakemake@wildcards$method
methods_config = fread(file=snakemake@params$methods_config, sep="\t")
nSteps = snakemake@params$nSteps %>% as.numeric()

outDir = snakemake@params$outDir
outFile = snakemake@output$outFile

methods_this = dplyr::filter(methods_config, method==method_this)
methods_this$binary = as.character(methods_this$binary)

if (methods_this$binary[1]=="TRUE"){
  print("binary")
  nSteps=2
} else {
  print("non-binary")
}

sampleKey_this = fread(methods_this$sampleKey[1], sep="\t")
sampleKey_this =  dplyr::select(sampleKey_this, biosample, GTExTissue) %>% dplyr::filter(GTExTissue!="")
sampleKey_this$key = paste0(sampleKey_this$GTExTissue, ".", sampleKey_this$biosample)

# first file
col_names = c("varChr", "varStart", "varEnd", "hgID", "GTExTissue", "eGene", "PIP", "TPM", "distance","enhChr", "enhStart", "enhEnd", "biosample", "TargetGene", "score" )
biosample_this = sampleKey_this$biosample[1]
GTExTissue_this = sampleKey_this$GTExTissue[1]
varIntFile = file.path(outDir, method_this, biosample_this, "GTExVariants-enhancerPredictionsInt.tsv.gz")
varInt = read.table(file=varIntFile, sep="\t", header=FALSE) %>% setNames(col_names)
varInt_filtered = dplyr::filter(varInt, GTExTissue==GTExTissue_this, eGene==TargetGene) %>%
  dplyr::select(varChr, varStart, varEnd, eGene, score) %>% distinct()
varInt_all = varInt_filtered

# the rest
if (nrow(sampleKey_this) > 1) {
  for (i in 2:length(sampleKey_this)){
    biosample_this=sampleKey_this$biosample[i]
    GTExTissue_this = sampleKey_this$GTExTissue[i]
    varIntFile = file.path(outDir, method_this, biosample_this, "GTExVariants-enhancerPredictionsInt.tsv.gz")
    varInt = read.table(file=varIntFile, sep="\t", header=FALSE) %>% setNames(col_names)
    varInt_filtered = dplyr::filter(varInt, GTExTissue==GTExTissue_this, eGene==TargetGene) %>%
      dplyr::select(varChr, varStart, varEnd, eGene, score) %>% distinct()
    varInt_all= rbind(varInt_all, varInt_filtered)
  }
}

scores = as.numeric(varInt_all$score)
prob_vector = seq(0, 1, length.out=nSteps)
thresholds = data.frame(unname(quantile(scores, na.rm=T, probs=prob_vector)))

colnames(thresholds) = "threshold"
write.table(thresholds, outFile, row.names=FALSE, col.names=FALSE, sep="\t")
