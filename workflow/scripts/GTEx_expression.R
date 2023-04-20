suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)})


main <- function() {
  expr.file = (snakemake@input$exprData) # highly expressed genes
  TPM.thresh = (snakemake@params$thresholdTPM)
  var.file = (snakemake@input$filteredGTExVar)
  out.file = (snakemake@output$GTExExpressedGenes)
  
  col.key = data.frame(GTExTissue=c('Brain_Cortex', 'Prostate', 'Muscle_Skeletal', 'Artery_Tibial', 
  'Skin_Not_Sun_Exposed_Suprapubic', 'Esophagus_Muscularis', 'Brain_Cerebellum', 'Pituitary', 
  'Adipose_Subcutaneous', 'Pancreas', 'Spleen', 'Adrenal_Gland', 'Testis', 'Lung', 
  'Skin_Sun_Exposed_Lower_leg', 'Heart_Atrial_Appendage', 'Adipose_Visceral_Omentum', 
  'Esophagus_Mucosa', 'Thyroid', 'Brain_Nucleus_accumbens_basal_ganglia', 'Colon_Sigmoid', 
  'Breast_Mammary_Tissue', 'Artery_Aorta', 'Cells_EBV-transformed_lymphocytes', 
  'Heart_Left_Ventricle', 'Artery_Coronary', 'Esophagus_Gastroesophageal_Junction', 
  'Nerve_Tibial', 'Liver', 'Whole_Blood', 'Colon_Transverse', 'Stomach'))
  col.key$colNum = c(15, 46, 41, 8, 47, 32, 14, 45, 3, 44, 50, 5, 52, 39, 48, 34, 4, 31, 53, 19, 28, 23, 6, 25, 35, 7, 30, 42, 38, 56, 29, 51)
  
  # read in expression data, select relevant columns, convert to long format with columns for tissue name and TPM
  expr = read.table(expr.file, header=FALSE, skip=3)
  expr = dplyr::select(expr, c(1, 2, col.key$colNum)) %>% setNames(c('eGene','hgnc',col.key$GTExTissue)) %>%
    pivot_longer(cols=-c(eGene, hgnc), names_to='GTExTissue', values_to='TPM')

  var = read.table(gzfile(var.file), header=FALSE) %>% setNames(c('chr','start','end','hgid','GTExTissue','eGene', 'PIP'))
  # merge variants and expr data, filter to TPM>threshold
  var.expr = left_join(var, expr, by=c('GTExTissue'='GTExTissue', 'eGene'='eGene')) %>% filter(TPM>TPM.thresh)
  var.expr$eGene = var.expr$hgnc
  var.expr = dplyr::select(var.expr, -hgnc)

  #out.file = file.path(outDir, "generalReference", "GTExVariants.filtered.PIP0.5.distalNoncoding.expressed.tsv")
  
  write.table(var.expr, file=out.file, sep='\t', row.names=FALSE, col.names=FALSE, quote=FALSE )
}


main()
