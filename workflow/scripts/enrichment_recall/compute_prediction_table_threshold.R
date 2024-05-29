suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(optparse)
  library(stringr)
  library(gdata)
  library(ggplot2)
  library(data.table)})

varIntFile = (snakemake@input$variantsPredictionsInt)
GTExVariantsFile = (snakemake@input$filteredGTExVariantsFinal)
thresholdSpan = read.table(snakemake@input$thresholdSpan, sep="\t", header=FALSE)
outFile = (snakemake@output$predTable)
GTExTissue.this = snakemake@wildcards$GTExTissue
Biosample.this = snakemake@wildcards$Biosample
method.this = snakemake@wildcards$method

varInt = read.table(file=varIntFile, sep="\t", header=FALSE) %>%
  setNames(c("varChr", "varStart", "varEnd", "variantID", "eGene", "GTExTissue", "PIP",  "distanceBin",
             "enhChr", "enhStart", "enhEnd", "Biosample", "TargetGene", "score"))
# filter to GTEx tissue, cell type (redundant)
varInt = dplyr::filter(varInt, GTExTissue==GTExTissue.this, Biosample==Biosample.this)

GTExVariants = read.table(file=GTExVariantsFile, sep="\t", header=FALSE) %>%
  setNames(c("varChr", "varStart", "varEnd", "variantID",  "eGene", "GTExTissue", "PIP", "distanceBin"))
GTExVariants = dplyr::filter(GTExVariants, GTExTissue==GTExTissue.this)
nVariantsTotal = dplyr::select(GTExVariants, varChr, varStart, varEnd, eGene) %>% distinct() %>% nrow()

# initialize pred table
predTable = data.frame(thresholdSpan) %>% setNames("threshold")
predTable$total.variants = nVariantsTotal
predTable$recall.total = 0 # fraction of variants overlapping enhancers
predTable$recall.linking = 0 # fraction of variants overlapping enhancers linked to correct gene
predTable$ correctGene.ifOverlap = 0 # fraction of variants linked to correct gene GIVEN they overlap enhancers 

# fill pred table
for (i in 1:nrow(predTable)){
  threshold.this = predTable$threshold[i]
  varInt.this = dplyr::filter(varInt, score>=threshold.this)

  nVariantsOverlappingEnhancers = dplyr::select(varInt.this, varChr, varStart, varEnd, eGene) %>% distinct() %>% nrow()
  nVariantsOverlappingEnhancersCorrectGene = dplyr::filter(varInt.this, eGene==TargetGene) %>% 
    dplyr::select(varChr, varStart, varEnd, eGene) %>% distinct() %>% nrow()
  
  predTable$total.variants[i] = nVariantsTotal
  predTable$recall.total[i] = nVariantsOverlappingEnhancers/nVariantsTotal
  predTable$recall.linking[i] = nVariantsOverlappingEnhancersCorrectGene/nVariantsTotal
  predTable$correctGene.ifOverlap[i] = nVariantsOverlappingEnhancersCorrectGene/nVariantsOverlappingEnhancers

}

predTable$Biosample = Biosample.this
predTable$GTExTissue = GTExTissue.this
predTable$method = method.this

write.table(predTable, file=outFile, sep="\t", quote=F, row.names=F, col.names=T)

