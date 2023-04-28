suppressPackageStartupMessages({
  library(biomaRt)
  library(dplyr)
  library(tidyr)
  library(optparse)
  library(stringr)
  library(gdata)
  library(ggplot2)})

varIntFile = (snakemake@input$variantsPredictionsInt)
GTExVariantsFile = (snakemake@input$filteredGTExVariantsFinal)
score.thresh = (snakemake@params$threshold) %>% as.numeric()
distance = (snakemake@params$distances) %>% as.character() %>% strsplit(" ") %>% unlist() %>% as.numeric()
outFile = (snakemake@output$predTable)
GTExTissue.this = snakemake@wildcards$GTExTissue
Biosample.this = snakemake@wildcards$Biosample


varInt = read.table(file=varIntFile, sep="\t", header=FALSE) %>%
  setNames(c("varChr", "varStart", "varEnd", "hgID", "GTExTissue", "eGene", "PIP", "TPM", "distance",
             "enhChr", "enhStart", "enhEnd", "Biosample", "TargetGene", "score"))
# filter to GTEx tissue, cell type (redundant), and score threshold
varInt = dplyr::filter(varInt, GTExTissue==GTExTissue.this, Biosample==Biosample.this, score>=score.thresh)
GTExVariants = read.table(file=GTExVariantsFile, sep="\t", header=FALSE) %>%
  setNames(c("varChr", "varStart", "varEnd", "hgID", "GTExTissue", "eGene", "PIP", "TPM", "distance"))
GTExVariants = dplyr::filter(GTExVariants, GTExTissue==GTExTissue.this)
  
# initialize pred table
predTable = data.frame(distance)
predTable$total.variants = 0
predTable$recall.total = 0 # fraction of variants overlapping enhancers
predTable$recall.linking = 0 # fraction of variants overlapping enhancers linked to correct gene
predTable$ correctGene.ifOverlap = 0 # fraction of variants linked to correct gene GIVEN they overlap enhancers 

# fill pred table
for (i in 1:nrow(predTable)){
  distance.this = predTable$distance[i]
  if (distance.this == 30000000){
    varInt.this = varInt
    GTExVariants.this = GTExVariants
  } else {
    varInt.this = dplyr::filter(varInt, distance==distance.this)
    GTExVariants.this = dplyr::filter(GTExVariants, distance==distance.this)
  }

  nVariantsTotal = dplyr::select(GTExVariants.this, varChr, varStart, varEnd, eGene) %>% distinct() %>% nrow()
  nVariantsOverlappingEnhancers = dplyr::select(varInt.this, varChr, varStart, varEnd, eGene) %>% distinct() %>% nrow()
  nVariantsOverlappingEnhancersCorrectGene = dplyr::filter(varInt.this, eGene==TargetGene) %>% 
    dplyr::select(varChr, varStart, varEnd, eGene) %>% distinct() %>% nrow()
  
  predTable$total.variants[i] = nVariantsTotal
  predTable$recall.total[i] = nVariantsOverlappingEnhancers/nVariantsTotal
  predTable$recall.linking[i] = nVariantsOverlappingEnhancersCorrectGene/nVariantsTotal
  predTable$correctGene.ifOverlap[i] = nVariantsOverlappingEnhancersCorrectGene/nVariantsOverlappingEnhancers

}
  
write.table(predTable, file=outFile, sep="\t", quote=F, row.names=F, col.names=T)

