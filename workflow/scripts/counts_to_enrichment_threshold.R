suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)})

## get files from snakemake
method = (snakemake@wildcards$method)
score.thresh = (snakemake@wildcards$threshold)
outDir = (snakemake@params$outDir)
countFile = (snakemake@input$countMatrix)
biosamples = (snakemake@params$biosamples) %>% (" ") %>% unlist() %>% sort()
varPerGTExTissueFile = (snakemake@input$variantsPerGTExTissue)
GTExTissues = (snakemake@params$GTExTissues) %>% strsplit(" ") %>% unlist %>% sort()
outFile = (snakemake@output$enrichmentTable)
totalCommonVar = (snakemake@params$nCommonVariants)

## read in files
# count matrix: each row = biosample, each col = tissue
countMatrix = read.table(file=countFile, header=TRUE, stringsAsFactors=FALSE)

# common variants by tissue/biosample
commonVarPerBiosample = data.frame(biosamples)%>% setNames(c("Biosample"))
commonVarPerBiosample$nCommonVariantsOverlappingEnhancers = 0
# read in per biosample and fill in df
# loop through biosamples
for (i in 1:nrow(biosamples)){
  sample.this = biosamples[i]
  commonVarPredIntFile = file.path(outDir, method, sample.this, "distalNoncodingBackgroundSNPs-enhancerPredictionsInt.tsv.gz")
  commonVarPredInt = read.table(file=commonVarPredIntFile, header=TRUE, stringsAsFactors=FALSE) %>%
    setNames(c("varChr", "varStart", "varEnd", "rsID", "enhChr", "enhStart", "enhEnd", "Biosample", "TargetGene", "score"))
  commonVarPredInt = dplyr::filter(commonVarPredInt, score>=score.thresh)
  
  counts.this = dplyr::select(commonVarPredInt, rsID) %>% distinct() %>% nrow()
  commonVarPerBiosample$nCommonVariantsOverlappingEnhancers[i] = counts.this
}

# variants per GTEx tissue
variantsByGTExTissue = read.table(varPerGTExTissueFile, header=TRUE, stringsAsFactors=FALSE) %>% 
  setNames("GTExTissue", "nVariantsGTExTissue")

# make enrichment matrix
enrMatrix = pivot_longer(countMatrix, cols=-Biosample, names_to='GTExTissue', values_to='nVariantsOverlappingEnhancers')
enrMatrix[enrMatrix=="Cells_EBV.transformed_lymphocytes"] = "Cells_EBV-transformed_lymphocytes"
enrMatrix = left_join(enrMatrix, variantsByGTExTissue, by='GTExTissue')
enrMatrix = left_join(enrMatrix, commonVarPerBiosample, by='Biosample')
enrMatrix$nCommonVariants = totalCommonVar
enrMatrix$enrichment = enrMatrix$nVariantsOverlappingEnhancers/enrMatrix$nVariantsGTExTissue/(enrMatrix$nCommonVariantsOverlappingEnhancers/enrMatrix$nCommonVariants)

write.table(enrMatrix, file=outFile, sep="\t", quote=F, row.names=F, col.names=T)
