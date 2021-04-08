suppressPackageStartupMessages({
  library(biomaRt)
  library(dplyr)
  library(tidyr)
  library(optparse)
  library(stringr)})

main <- function() {
  option_list <- list(
    make_option(c("--genes"), type="character", default=NA, help="file of genes (HGNC symbol) to filter input to; cols = chr,start,end,gene"),
    make_option(c("--col"), type="numeric", default=4, help="col # of input file containing ensembl IDs"))
  opt = parse_args(OptionParser(option_list=option_list))
  
  col = opt$col; ABC.gene.file=opt$genes

    # read ABC data
  ABC.genes = read.csv(ABC.gene.file, sep='\t', header=FALSE)
  if (ncol(ABC.genes)>1){
    ABC.genes = dplyr::select(ABC.genes, chr='V1',start='V2',end='V3',hgnc.IDs='V4')
  } else {
    colnames(ABC.genes) = 'hgnc.IDs'
  }
  # read dataframe, remove values after decimal in ens ID column
  df = readLines(con = file("stdin")) %>% as.data.frame()

  nCols = str_count(df[1,], '\t') + 1
  df = separate(df, col='.',into=as.character(1:nCols),sep="\t")
  
  df[[col]] = substr(df[[col]], 1,15) # get rid of decimals after ensembl IDs
  
  # convert ABC genes to ens id
  ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
  ens.IDs = getBM(attributes=c('ensembl_gene_id','hgnc_symbol'), filters = 'hgnc_symbol', values = ABC.genes$hgnc.IDs, mart = ensembl)
  ABC.genes = left_join(ABC.genes, ens.IDs, by=c("hgnc.IDs"="hgnc_symbol"))
  
  # filter dataset to only having IDs in ABC list
  colnames(df)[col] = 'gene.IDs'
  df = dplyr::filter(df, is.element(df$gene.IDs, ABC.genes$ensembl_gene_id))
  
  # write table
  
  write.table(df, file="", sep="\t", quote=F, row.names=F, col.names=F)}

main()
