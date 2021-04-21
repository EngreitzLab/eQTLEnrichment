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
    make_option(c("--col"), type="numeric", default=4, help="col # of input file containing ensembl IDs"),
    make_option(c("--id"), type="character", default="ensembl", help="type of gene ids in file being filtered ('hgnc' or 'ensembl'), defaults to 'ensembl'"))
  
  opt = parse_args(OptionParser(option_list=option_list))
  inFile=opt$input; col = opt$col; ABC.gene.file=opt$genes; idType=opt$id
  
  # read dataframe
  #df = readLines(con = file("stdin")) %>% as.data.frame()
  #nCols = str_count(df[1,], '\t') + 1
  #df = separate(df, col='.',into=as.character(1:nCols),sep="\t")

  df = read.table(inFile, header=FALSE)
  
  # read ABC data
  ABC.genes = read.table(ABC.gene.file, header=FALSE)
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
  
  # write table
  write.table(df, file="", sep="\t", quote=F, row.names=F, col.names=F)} else {
    
    colnames(df)[col] = 'gene.IDs'
    df = dplyr::filter(df, is.element(df$gene.IDs, ABC.genes$hgnc.IDs))
    write.table(df, file="", sep="\t", quote=F, row.names=F, col.names=F)}
    
  }
  
  


main()
