suppressPackageStartupMessages({
  library(biomaRt)
  library(dplyr)
  library(tidyr)
  library(optparse)
  library(stringr)})


main <- function() {
  option_list <- list(
    make_option(c("--variants"), type="character", default=NA, help="tsv with variants"),
    make_option(c("--pred"), type="character", default=NA, help="intersection of variants with enhancer predictions (col 1 = unique ID, col 2 = predicted target gene"),
    make_option(c("--genes"), type="character", default=NA, help="list of expressed genes for this tissue/biosample"))
  
  
  opt = parse_args(OptionParser(option_list=option_list))
  varFile = opt$variants; predictionsInt = opt$pred; exprGenes = opt$genes
  
  variants = read.table(file=varFile,header=TRUE, stringsAsFactor=FALSE); 
  predictionsInt = read.table(file=predictionsInt,header=FALSE, stringsAsFactor=FALSE)
  colnames(predictionsInt) = c('unique_id','targetGene')
  exprGenes = read.table(file=exprGenes, header=FALSE, stringsAsFactor=FALSE)
  colnames(exprGenes) = c('chr','start','end','gene','n','dir')

    
  for (i in 1:nrow(variants)) {
    temp = filter(predictionsInt, unique_id==variants$unique_id[i]) %>% 
      filter(is.element(targetGene, exprGenes$gene))
    
    if (nrow(temp)==0){
      variants$predictionClass[i] = 'noOverlap'
      variants$nEnhancers[i]=0
    } else if (is.element(variants$eGene[i], temp$targetGene)) {
      variants$predictionClass[i] = 'inEnhancer-correctGene'
      variants$nEnhancers[i] = nrow(temp)
    } else {
      variants$predictionClass[i] = 'inEnhancer-incorrectGene'
      variants$nEnhancers[i] = nrow(temp)
    }
  }
  
  # filter out rows with eGenes connected to un-expressed genes
  variants = filter(variants, is.element(eGene, exprGenes$gene))
  write.table(variants, file="", sep="\t", quote=F, row.names=F, col.names=T)
}

main()
