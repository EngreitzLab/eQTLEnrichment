suppressPackageStartupMessages({
  library(biomaRt)
  library(dplyr)
  library(tidyr)
  library(optparse)
  library(stringr)
  library(gdata)})


main <- function() {
  option_list <- list(
    make_option(c("--variants"), type="character", default=NA, help="tsv with variants"),
    make_option(c("--TSSint"), type="character", default=NA, help="intersection of variants with TSS (col 1 = unique ID, col 2 = gene with nearby TSS"))
  
  opt = parse_args(OptionParser(option_list=option_list))
  varFile = opt$variants; intFile = opt$TSSint
  
  variants = read.csv(file=varFile, sep='\t',header=FALSE, stringsAsFactor=FALSE, fill=TRUE); 
  colnames(variants)=c('chr','start','end','unique_id','eGene','PIP','closestGene','closestTSS')
  TSSint = read.csv(file=intFile, sep='\t',header=FALSE, stringsAsFactor=FALSE, fill=TRUE)
  colnames(TSSint) = c('unique_id','closeTSS')
  
  for (i in 1:nrow(variants)) {
    temp = filter(TSSint, unique_id==variants$unique_id[i])
    #variants$nearbyTSS[i] = is.element(variants$eGene[i], temp$closeTSS)
    variants$nearbyTSS[i] = (any(gdata::startsWith(str=temp$closeTSS, pattern=variants$eGene[i])))
    variants$nTSSnear[i] = nrow(temp)
  }
  
  write.table(variants, file="", sep="\t", quote=F, row.names=F, col.names=T)
}

main()
