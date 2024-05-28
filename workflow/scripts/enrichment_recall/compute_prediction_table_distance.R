suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(gdata)})

varIntFile = (snakemake@input$variantsPredictionsInt)
GTExVariantsFile = (snakemake@input$filteredGTExVariantsFinal)
score.thresh = (snakemake@params$threshold) %>% as.numeric()
distances_min = (snakemake@params$distances_min) %>% as.character() %>% strsplit(" ") %>% unlist() %>% as.numeric()
distances_max = (snakemake@params$distances_max) %>% as.character() %>% strsplit(" ") %>% unlist() %>% as.numeric()
outFile = (snakemake@output$predTable)
GTExTissue.this = snakemake@wildcards$GTExTissue
Biosample.this = snakemake@wildcards$Biosample
method.this = snakemake@wildcards$method


varInt = read.table(file=varIntFile, sep="\t", header=FALSE) %>%
  setNames(c("varChr", "varStart", "varEnd", "variantID", "eGene", "GTExTissue", "PIP", "distanceBin",
             "enhChr", "enhStart", "enhEnd", "Biosample", "TargetGene", "score"))
# filter to GTEx tissue, cell type (redundant), and score threshold
varInt = dplyr::filter(varInt, GTExTissue==GTExTissue.this, Biosample==Biosample.this, score>=score.thresh)
GTExVariants = read.table(file=GTExVariantsFile, sep="\t", header=FALSE) %>%
  setNames(c("varChr", "varStart", "varEnd", "variantID", "eGene", "GTExTissue", "PIP", "distanceBin"))
GTExVariants = dplyr::filter(GTExVariants, GTExTissue==GTExTissue.this)
  
# initialize pred table

predTable = data.frame(distance_min=distances_min, distance_max=distances_max)
predTable$total.variants = 0
predTable$recall.total = 0 # fraction of variants overlapping enhancers
predTable$recall.linking = 0 # fraction of variants overlapping enhancers linked to correct gene
predTable$ correctGene.ifOverlap = 0 # fraction of variants linked to correct gene GIVEN they overlap enhancers 

# fill pred table
for (i in 1:nrow(predTable)){
  distance.this = predTable$distance_max[i]
  if (distance.this== 30000){ # all variants
    varInt.this = varInt
    GTExVariants.this = GTExVariants
  } else {
    varInt.this = dplyr::filter(varInt, distanceBin==distance.this)
    GTExVariants.this = dplyr::filter(GTExVariants, distanceBin==distance.this)
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

predTable$Biosample = Biosample.this
predTable$GTExTissue = GTExTissue.this
predTable$method = method.this
  
write.table(predTable, file=outFile, sep="\t", quote=F, row.names=F, col.names=T)

