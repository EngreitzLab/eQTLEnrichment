suppressPackageStartupMessages({
  library(biomaRt)
  library(dplyr)
  library(tidyr)
  library(optparse)
  library(stringr)
  library(gdata)
  library(ggplot2)})

  option_list <- list(
    make_option(c("--variants"), type="character", default=NA, help="tsv with variants"),
    make_option(c("--pred"), type="character", default=NA, help="intersection of variants with enhancer predictions (col 1 = unique ID, col 2 = predicted target gene, col 3 = score"),
    make_option(c("--outDir"), type="character", default=NA, help="output directory for PR plots"),
    make_option(c("--outThresh"), type="character", default=NA, help="outut file for threshold table"),
    make_option(c("--biosample"), type="character", default=NA, help="biosample for enhancer predictions"))
  
  opt = parse_args(OptionParser(option_list=option_list))
  varFile = opt$variants; predictionsIntFile = opt$pred; outDir = opt$outDir; biosample=opt$biosample
  
  variants = read.table(file=varFile,header=FALSE, stringsAsFactor=FALSE); 
  colnames(variants) = c("chr", "start", "end", "unique_id", "GTExTissue", "eGene", "PIP", "TPM")
  #variants$annotations=""
  
  predIntSize = file.info(predictionsIntFile)$size
  size.threshold=50
  if (predIntSize<size.threshold){
    predictionsInt=matrix(nrow=1,ncol=3) %>% data.frame()
  } else {
    predictionsInt = read.table(file=predictionsIntFile,header=FALSE, stringsAsFactor=FALSE)
  }
  colnames(predictionsInt) = c('unique_id', 'targetGene', 'Score')
  predictionsInt$Score=as.numeric(predictionsInt$Score)
  
  tissue = str_split(predictionsIntFile, '/')[[1]]; tissue = tissue[length(tissue)] %>% str_split(pattern='\\.'); tissue = tissue[[1]][1]
  method = str_split(predictionsIntFile, '/')[[1]]; method = method[length(method)-2]

write_table <- function() {    
  for (i in 1:nrow(variants)) {
    temp = filter(predictionsInt, unique_id==variants$unique_id[i])
    variants$nEnhancers[i] = nrow(temp)
    
    if (nrow(temp)==0 && !(variants$eGene[i] %in% temp$targetGene)) {
      variants$predictionClass[i] = 'noOverlap'
      if  (method=="reads_by_dist_to_tss" || method=="reads_by_dist_to_tss_norm") {
        variants$Score[i] = 10
      } else if (method=="dist_to_gene" || method=="dist_to_tss") {
        variants$Score[i] = 0
      } else {
        variants$Score[i] = 1
      }
    }
    
    if (nrow(temp)>0 && (variants$eGene[i] %in% temp$targetGene)) {
      variants$predictionClass[i] = 'inEnhancer-correctGene'
      if  (method=="reads_by_dist_to_tss" || method=="reads_by_dist_to_tss_norm") {
        variants$Score[i]=max(temp$Score, na.rm=TRUE)
      }  else if (method=="dist_to_gene" || method=="dist_to_tss") {
        variants$Score[i]=min(temp$Score)
      } else {
        variants$Score[i]=1
      }
    } 
    
    if (nrow(temp)>0 && !(variants$eGene[i] %in% temp$targetGene)) {
      variants$predictionClass[i] = 'inEnhancer-incorrectGene'
      variants$Score[i]=if_else(length(min(temp$Score))==1, min(temp$Score), min(temp$Score)[1])
    }
  }
  
  write.table(variants, file="", sep="\t", quote=F, row.names=F, col.names=T)
}


write_table()
