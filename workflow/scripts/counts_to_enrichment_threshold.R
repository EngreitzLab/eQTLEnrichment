suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)})

## get files from snakemake
method = (snakemake@wildcards$method)
score.thresh = (snakemake@wildcards$threshold) %>% as.numeric()
print(score.thresh)
outDir = (snakemake@params$outDir)
countFile = (snakemake@input$countMatrix)
biosamples = (snakemake@params$biosamples) %>% strsplit(" ") %>% unlist() %>% sort()
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
for (i in 1:length(biosamples)){
  sample.this = biosamples[i]
  print(sample.this)
  commonVarPredIntFile = file.path(outDir, method, sample.this, "distalNoncodingBackgroundSNPs-enhancerPredictionsInt.tsv.gz")
  commonVarPredInt = read.table(file=commonVarPredIntFile, header=TRUE, stringsAsFactors=FALSE) %>%
    setNames(c("varChr", "varStart", "varEnd", "rsID", "enhChr", "enhStart", "enhEnd", "Biosample", "TargetGene", "score"))
  print(nrow(commonVarPredInt))
  commonVarPredInt = dplyr::filter(commonVarPredInt, score>=score.thresh)
  print(nrow(commonVarPredInt))
  
  counts.this.df = dplyr::select(commonVarPredInt, rsID) %>% unique()
  counts.this = nrow(counts.this.df)
  print(tail(counts.this.df))
  print(counts.this)
  commonVarPerBiosample$nCommonVariantsOverlappingEnhancers[i] = counts.this
}

# variants per GTEx tissue
variantsByGTExTissue = read.table(varPerGTExTissueFile, header=TRUE, stringsAsFactors=FALSE) %>% 
  setNames(c("GTExTissue", "nVariantsGTExTissue"))

# make enrichment matrix
enrMatrix = pivot_longer(countMatrix, cols=-Biosample, names_to='GTExTissue', values_to='nVariantsOverlappingEnhancers')
enrMatrix[enrMatrix=="Cells_EBV.transformed_lymphocytes"] = "Cells_EBV-transformed_lymphocytes"
enrMatrix = left_join(enrMatrix, variantsByGTExTissue, by='GTExTissue')
enrMatrix = left_join(enrMatrix, commonVarPerBiosample, by='Biosample')
enrMatrix$nCommonVariants = totalCommonVar
enrMatrix$enrichment = enrMatrix$nVariantsOverlappingEnhancers/enrMatrix$nVariantsGTExTissue/(enrMatrix$nCommonVariantsOverlappingEnhancers/enrMatrix$nCommonVariants)

write.table(enrMatrix, file=outFile, sep="\t", quote=F, row.names=F, col.names=T)
