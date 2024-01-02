suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)})

## get files from snakemake
method = (snakemake@wildcards$method)
thresholdSpan = (read.table(snakemake@input$thresholdSpan, sep='\t', header=FALSE)) %>% setNames("threshold")
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
commonVarPerBiosample = data.frame(biosamples)%>% setNames(c("Biosample"))
commonVarPerBiosample$nCommonVariantsOverlappingEnhancers = 0

for (k in 1:nrow(thresholdSpan)){
  score.thresh=thresholdSpan$threshold[k]

  # common variants by tissue/biosample
  # read in per biosample and fill in df
  # loop through biosamples
  for (i in 1:length(biosamples)){
    sample.this = biosamples[i]
    commonVarPredIntFile = file.path(outDir, method, sample.this, "distalNoncodingBackgroundSNPs-enhancerPredictionsInt.tsv.gz")
    commonVarPredInt = read.table(file=commonVarPredIntFile, header=TRUE, stringsAsFactors=FALSE) %>%
      setNames(c("varChr", "varStart", "varEnd", "rsID", "enhChr", "enhStart", "enhEnd", "Biosample", "TargetGene", "score"))
    commonVarPredInt = dplyr::filter(commonVarPredInt, score>=score.thresh)
    
    counts.this.df = dplyr::select(commonVarPredInt, rsID) %>% unique()
    counts.this = nrow(counts.this.df)
    commonVarPerBiosample$nCommonVariantsOverlappingEnhancers[i] = counts.this
  }
  # variants per GTEx tissue
  variantsByGTExTissue = read.table(varPerGTExTissueFile, header=TRUE, stringsAsFactors=FALSE) %>% 
    setNames(c("GTExTissue", "nVariantsGTExTissue"))
  # make enrichment matrix
  enrMatrix = pivot_longer(countMatrix, cols=-c(Biosample,threshold), names_to='GTExTissue', values_to='nVariantsOverlappingEnhancers')
  enrMatrix = dplyr::filter(enrMatrix, threshold==score.thresh)
  enrMatrix[enrMatrix=="Cells_EBV.transformed_lymphocytes"] = "Cells_EBV-transformed_lymphocytes"
  enrMatrix = left_join(enrMatrix, variantsByGTExTissue, by='GTExTissue')
  enrMatrix = left_join(enrMatrix, commonVarPerBiosample, by='Biosample')
  enrMatrix$nCommonVariants = totalCommonVar
  enrMatrix$enrichment = enrMatrix$nVariantsOverlappingEnhancers/enrMatrix$nVariantsGTExTissue/(enrMatrix$nCommonVariantsOverlappingEnhancers/enrMatrix$nCommonVariants)
  
  enrMatrix$threshold = score.thresh
  if (k==1){
    enrMatrix_all = enrMatrix
  } else {
    enrMatrix_all = rbind(enrMatrix_all, enrMatrix)
  }
}

## stats about risk ratio (RR) aka enrichment
# calculate CI of RR and SE(log RR); see: https://sphweb.bumc.bu.edu/otlt/mph-modules/bs/bs704_confidence_intervals/bs704_confidence_intervals8.html
z = 1.96 # for 95% CI
calcs = enrMatrix_all
calcs$n1 = calcs$nVariantsGTExTissue
calcs$x1 = calcs$nVariantsOverlappingEnhancers
calcs$n2 = calcs$nCommonVariants
calcs$x2 = calcs$nCommonVariantsOverlappingEnhancers

calcs$SE_log_enr = with(calcs, sqrt(((n1-x1)/x1)/n1) + ((n2-x2)/x2)/n2)
calcs$log_enr = log(calcs$enrichment)
calcs$log_CI_enr_low = with(calcs, log_enr - z*SE_log_enr)
calcs$log_CI_enr_high = with(calcs, log_enr + z*SE_log_enr)
calcs$CI_enr_low = exp(calcs$log_CI_enr_low)
calcs$CI_enr_high = exp(calcs$log_CI_enr_high)

enrMatrix_all$CI_enr_low = calcs$CI_enr_low
enrMatrix_all$CI_enr_high = calcs$CI_enr_high
enrMatrix_all$SE_log_enr = calcs$SE_log_enr
  
## write output
write.table(enrMatrix_all, file=outFile, sep="\t", quote=F, row.names=F, col.names=T)
