suppressPackageStartupMessages({
  library(biomaRt)
  library(dplyr)
  library(tidyr)
  library(optparse)
  library(stringr)})

main <- function() {
  option_list <- list(
    make_option(c("--input"), type="character", default=NA, help="input file"),
    make_option(c("--genes"), type="character", default=NA, help="file of genes (HGNC symbol) to filter input to; cols = chr,start,end,gene"),
    make_option(c("--id_col"), type="numeric", default=4, help="col # of input file containing ensembl IDs"),
    make_option(c("--id"), type="character", default="ensembl", help="type of gene ids in file being filtered ('hgnc' or 'ensembl'), defaults to 'ensembl'"),
    make_option(c("--biosample"), type="character", default = NA, help="biosample for this prediction file"),
    make_option(c("--biosample_col"), type="numeric", default = 4, help="col # with biosample"),
    make_option(c("--invert"), type="character", default = "False", help="invert score?"),
    make_option(c("--score_col"), type="numeric", default = 6, help="col # with predictor score"))
  
  opt = parse_args(OptionParser(option_list=option_list))
  inFile=opt$input; col = opt$id_col; ABC.gene.file=opt$genes; idType=opt$id
  biosample=opt$biosample; biosample.col = opt$biosample_col
  invert=opt$invert; score.col=opt$score_col

  df = read.table(inFile, header=FALSE, fill=TRUE)

  # read ABC data
  ABC.genes = read.table(ABC.gene.file, header=FALSE, fill=TRUE)
  if (ncol(ABC.genes)>1){
    ABC.genes = dplyr::select(ABC.genes, chr='V1',start='V2',end='V3',hgnc.IDs='V4')
  } else {
    colnames(ABC.genes)[1] = 'hgnc.IDs'
  }

  if(idType=='ensembl'){
    df[[col]] = substr(df[[col]], 1,15) # get rid of decimals after ensembl IDs
  
    # convert ABC genes to ens id
    ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
    ens.IDs = getBM(attributes=c('ensembl_gene_id','hgnc_symbol'), filters = 'hgnc_symbol', values = ABC.genes$hgnc.IDs, mart = ensembl)
    ABC.genes = left_join(ABC.genes, ens.IDs, by=c("hgnc.IDs"="hgnc_symbol"))
    
    # filter dataset to only having IDs in ABC list
    colnames(df)[col] = 'gene.IDs'
    df = dplyr::filter(df, is.element(df$gene.IDs, ABC.genes$ensembl_gene_id))
  
} else {
    colnames(df)[col] = 'gene.IDs'
    df = dplyr::filter(df, is.element(df$gene.IDs, ABC.genes$hgnc.IDs))
  }
  
  # set biosample name if necessary
  if (!is.na(biosample)){
    colnames(df)[biosample.col] = 'biosample'
    df$biosample = biosample
  }
  
  # invert score if necessary
  if (invert == "True"){
    colnames(df)[score.col] = "score"
    df$score = -df$score
  }
  
  # write table
  write.table(df, file="", sep="\t", quote=F, row.names=F, col.names=F)
  }
  
  


main()
