suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(data.table)})

## get files from snakemake
method = (snakemake@wildcards$method)
distance_this = (snakemake@wildcards$distance_max) %>% as.numeric()
distances_max = (snakemake@params$distances_max) %>% as.character() %>% strsplit(" ") %>% as.numeric()
score.thresh = (snakemake@params$threshold) %>% as.numeric()
countFile = (snakemake@input$countMatrix)
bgCommonVar_file = (snakemake@input$commonVarCount)
commonVarInt_files = (snakemake@input$commonVarInt)
biosamples = (snakemake@params$biosamples) %>% strsplit(" ") %>% unlist()
varPerGTExTissueFile = (snakemake@input$variantsPerGTExTissueByDist)
outFile = (snakemake@output$enrichmentTable)

# count matrix: each row = biosample, each col = tissue
countMatrix = read.table(file=countFile, header=TRUE, stringsAsFactors=FALSE)

# common variants by tissue/biosample
commonVarPerBiosample = data.frame(biosamples)%>% setNames(c("Biosample"))
commonVarPerBiosample$nCommonVariantsOverlappingEnhancers = 0

for (i in 1:length(biosamples)){
  sample.this = biosamples[i]
  commonVarPredIntFile = commonVarInt_files[i]
  commonVarPredInt = read.table(file=commonVarPredIntFile, header=TRUE, stringsAsFactors=FALSE, sep="\t") %>%
    setNames(c("varChr", "varStart", "varEnd", "rsID", "enhChr", "enhStart", "enhEnd", "Biosample", "TargetGene", "score"))
  counts.this = dplyr::select(commonVarPredInt, rsID) %>% unique() %>% nrow()

  commonVarPerBiosample$nCommonVariantsOverlappingEnhancers[i] = counts.this
}

# variants per GTEx tissue
variantsByGTExTissue_init = read.table(varPerGTExTissueFile, header=TRUE, stringsAsFactors=FALSE, sep="\t") # col names: tissue, n_all, n_bin1, etc.
col_name = ifelse(distance_this==30000, "n_all", paste0("n_bin", which(distances_max==distance_this)))
variantsByGTExTissue = data.frame(GTExTissue = variantsByGTExTissue_init[["tissue"]], nVariantsGTExTissue = variantsByGTExTissue_init[[col_name]])

# background common variants
bgCommonVar_file = (snakemake@input$commonVarCount)
bgCommonVar = read.table(bgCommonVar_file, header=FALSE, sep="\t") %>% setNames("n")
bgCommonVar = as.numeric(bgCommonVar$n)

# make enrichment matrix
enrMatrix = pivot_longer(countMatrix, cols=-Biosample, names_to='GTExTissue', values_to='nVariantsOverlappingEnhancers')
enrMatrix[enrMatrix=="Cells_EBV.transformed_lymphocytes"] = "Cells_EBV-transformed_lymphocytes"
enrMatrix = left_join(enrMatrix, variantsByGTExTissue, by='GTExTissue')
enrMatrix = left_join(enrMatrix, commonVarPerBiosample, by='Biosample')
enrMatrix$nCommonVariants = bgCommonVar
enrMatrix$enrichment = enrMatrix$nVariantsOverlappingEnhancers/enrMatrix$nVariantsGTExTissue/(enrMatrix$nCommonVariantsOverlappingEnhancers/enrMatrix$nCommonVariants)

## stats about risk ratio (RR) aka enrichment
# calculate CI of RR and SE(log RR); see: https://sphweb.bumc.bu.edu/otlt/mph-modules/bs/bs704_confidence_intervals/bs704_confidence_intervals8.html
z = 1.96 # for 95% CI
calcs = enrMatrix
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

enrMatrix$CI_enr_low = calcs$CI_enr_low
enrMatrix$CI_enr_high = calcs$CI_enr_high
enrMatrix$SE_log_enr = calcs$SE_log_enr

enrMatrix$method = method
enrMatrix$distance_min = snakemake@wildcards$distance_min
enrMatrix$distance_max = snakemake@wildcards$distance_max

write.table(enrMatrix, file=outFile, sep="\t", quote=F, row.names=F, col.names=T)
