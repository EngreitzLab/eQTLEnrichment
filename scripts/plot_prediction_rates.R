suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(optparse)
  library(stringr)
  library(ggplot2)})

main <- function() {
  option_list <- list(
    make_option(c("--tables"), type="character", default=NA, help="tsv with variants"),
    make_option(c("--out"), type="character", default=NA, help="file path for plot"))
    opt = parse_args(OptionParser(option_list=option_list))
    outFile = opt$out
    tables = strsplit(opt$tables, " ") %>% unlist()
    
    # get list of GTEx tissues represented
    GTExTissues = str_split(tables, '/') %>% data.frame() %>% t()
    GTExTissues = GTExTissues[,ncol(GTExTissues)] %>% str_split('\\.') %>% data.frame() %>% t() %>% data.frame() %>% dplyr::select(X1) %>% distinct()

    # add rates for proximity-based E-G predictions
    df = data.frame(matrix(ncol=4,nrow=0, dimnames=list(NULL)))
    colnames(df) = c("GTExTissue", "metric", "method","value")
    
    for (x in GTExTissues$X1) {
      # select a table
      
      temp = tables[str_detect(tables, x)]
      temp = read.table(file=temp[1],header=TRUE,stringsAsFactors=FALSE)
      
      # calculate rates
      closest.gene.rate = nrow(filter(temp, closestGene==eGene))/nrow(temp)
      closest.TSS.rate = nrow(filter(temp, closestTSS==eGene))/nrow(temp)
      nearby.TSS.rate = nrow(filter(temp, nearbyTSS==TRUE))/nrow(temp)
      
      # add to df
      df = add_row(df, GTExTissue=x, metric="Closest gene body", method="Proximity", value=closest.gene.rate) %>%
        add_row(GTExTissue=x, metric="Closest TSS", method="Proximity", value=closest.TSS.rate) %>%
        add_row(GTExTissue=x, metric="TSS within 100 kb", method="Proximity", value=nearby.TSS.rate)
    }
    colnames(df) = c("GTExTissue", "metric", "method","value")

    # read in prediction rates from each method
    for (fileName in tables){
      predTable = read.table(file=fileName, header=TRUE, stringsAsFactors=FALSE)
      tissue = str_split(fileName, '/')[[1]]; tissue = tissue[length(tissue)] %>% str_split(pattern='\\.'); tissue = tissue[[1]][1]
      method = str_split(fileName, '/')[[1]]; method = method[length(method)-1]
      
      prediction.rate.GivenEnhancer = nrow(filter(predTable, predictionClass=='inEnhancer-correctGene'))/nrow(filter(predTable,predictionClass!='noOverlap'))
      prediction.rate.total = nrow(filter(predTable, predictionClass=='inEnhancer-correctGene'))/nrow(predTable)
      prediction.rate.inEnhancer=nrow(filter(predTable,predictionClass!='noOverlap'))/nrow(predTable)
      
      df = rbind(df, c(tissue, "Prediction rate given\nvariant in enhancer", method, prediction.rate.GivenEnhancer)) %>%
        rbind(c(tissue, "% variants in\nany enhancer", method, prediction.rate.inEnhancer)) %>%
        rbind(c(tissue, "Prediction rate", method, prediction.rate.total))
    }
    
    df$value=as.numeric(df$value)
    df$GTExTissue=as.factor(df$GTExTissue)
    df$metric = factor(df$metric); df$metric=factor(df$metric, levels = levels(df$metric)[c(2,3,6,1,4,5)])

    # graph
    g = ggplot(df, aes(x=metric, y=value, fill=method)) + geom_bar(stat="identity", position="dodge", width=0.5) + facet_grid(GTExTissue~.) + theme_minimal() + theme(axis.text.x = element_text(angle = 60,hjust=1)) + ylim(0,1)
    
    pdf(file=outFile,width=8,height=nrow(GTExTissues)*2); print(g); dev.off()

}

main()
