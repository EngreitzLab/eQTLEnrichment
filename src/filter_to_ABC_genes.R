library(biomaRt)
library(dplyr)
library(tidyr)
library(optparse)

main <- function() {
  option_list <- list(
    make_option(c("--genes"), type="character", default=NA, help="file of genes (HGNC symbol) to filter input to; cols = chr,start,end,gene"),
    make_option(c("--col"), type="numeric", default=4, help="col # of input file containing ensembl IDs"))
  opt = parse_args(OptionParser(option_list=option_list))
  
  colID = opt$col; ABC.gene.file=opt$genes
  
  # read ABC data
  ABC.genes = read.csv(ABC.gene.file, sep='\t', header=FALSE) %>% dplyr::select(chr='V1',start='V2',end='V3',hgnc.IDs='V4')
  
  # read dataframe, remove values after decimal in ens ID column
  df = readLines(con = file("stdin")) %>% as.data.frame() %>% separate(col='.',into=c('chr','start','end','hgid','tissue','ens.IDs','PIP'),sep='\t')
  df$ens.IDs = substr(df$ens.IDs, 1,15)
  
  # convert ABC genes to ens id
  ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
  ens.IDs = getBM(attributes=c('ensembl_gene_id','hgnc_symbol'), filters = 'hgnc_symbol', values = ABC.genes$hgnc.IDs, mart = ensembl)
  ABC.genes = left_join(ABC.genes, ens.IDs, by=c("hgnc.IDs"="hgnc_symbol"))
  
  # filter dataset to only having IDs in ABC list
  df = dplyr::filter(df, is.element(ens.IDs, ABC.genes$ensembl_gene_id))
  
  # write table
  
  write.table(df, file="", sep="\t", quote=F, row.names=F, col.names=F)}

main()
