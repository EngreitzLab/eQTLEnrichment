suppressPackageStartupMessages({
  library(biomaRt)
  library(dplyr)
  library(tidyr)
  library(optparse)
  library(stringr)})


main <- function() {
  option_list <- list( 
    make_option(c("--col"), type="numeric", default=5, help="col # of input file (read by stdin) containing ensembl IDs"))
    opt = parse_args(OptionParser(option_list=option_list))
    col = opt$col

    df = readLines(con = file("stdin")) %>% as.data.frame()
    nCols = str_count(df[1,], '\t') + 1
    df = separate(df, col='.',into=as.character(1:nCols),sep="\t")

    df[col] = substr(df[[col]], 1,15) # get rid of decimals after ensembl IDs

    ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
    key = getBM(attributes=c('ensembl_gene_id','hgnc_symbol'), filters = 'ensembl_gene_id', values = df[col], mart = ensembl)
    colnames(key) = c(as.character(col), 'hgnc_symbol')
    df = left_join(df, key)

    df[col] = df$hgnc_symbol
    df = filter(df, hgnc_symbol!="") # eliminate rows with no matching hgnc symbol
    df = dplyr::select(df, -hgnc_symbol)

    write.table(df, file="", sep="\t", quote=F, row.names=F, col.names=F)
    }

main()
