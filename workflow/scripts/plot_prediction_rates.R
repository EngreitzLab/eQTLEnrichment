suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(optparse)
  library(stringr)
  library(ggplot2)
  library(viridis)
  library(plotly)
  library(plyr)})

main <- function() {
  option_list <- list(
    make_option(c("--tables"), type="character", default=NA, help="tsv with variants"),
    make_option(c("--out"), type="character", default=NA, help="output directory"),
    make_option(c("--colors"), type="character", default=NA, help=".rds for color palette"))
  
    opt = parse_args(OptionParser(option_list=option_list))
    outDir = opt$out; outFile.pred = paste0(outDir, "/predictionRates.pdf"); outFile.PPV = paste0(outDir, "/PPV.pdf")
    tables = strsplit(opt$tables, " ") %>% unlist()
    cpList = readRDS(opt$colors)

    # get list of GTEx tissues represented
    GTExTissues = str_split(tables, '/') %>% data.frame() %>% t()
    GTExTissues = GTExTissues[,ncol(GTExTissues)] %>% str_split('\\.') %>% data.frame() %>% t() %>% data.frame() %>% dplyr::select(X1) %>% distinct()

    # add rates for proximity-based E-G predictions
    df = data.frame('temp','temp','temp',1)
    colnames(df) = c("GTExTissue", "metric", "method","value")

    
    # read in prediction rates from each method
    for (fileName in tables){
      predTable = read.table(file=fileName, header=TRUE, stringsAsFactors=FALSE)
      tissue = str_split(fileName, '/')[[1]]; tissue = tissue[length(tissue)] %>% str_split(pattern='\\.'); tissue = tissue[[1]][1]
      method = str_split(fileName, '/')[[1]]; method = method[length(method)-2]
      
      prediction.rate.GivenEnhancer = nrow(filter(predTable, predictionClass=='inEnhancer-correctGene'))/nrow(filter(predTable,predictionClass!='noOverlap'))
      prediction.rate.total = nrow(filter(predTable, predictionClass=='inEnhancer-correctGene'))/nrow(predTable)
      prediction.rate.inEnhancer=nrow(filter(predTable,predictionClass!='noOverlap'))/nrow(predTable)
      denom = dplyr::select(predTable, chr, start, end, nEnhancers) %>% distinct() 
      PPV = nrow(filter(predTable, predictionClass=='inEnhancer-correctGene'))/(sum(denom$nEnhancers))
      
      df = rbind(df, c(tissue, "Prediction rate given \nvariant in enhancer", method, prediction.rate.GivenEnhancer)) %>%
        rbind(c(tissue, "% variants in \nany enhancer", method, prediction.rate.inEnhancer)) %>%
        rbind(c(tissue, "Prediction rate", method, prediction.rate.total)) %>%
        rbind(c(tissue, "Positive predictive value", method, PPV))
    }
    
    df = filter(df, GTExTissue!="temp")
    
    df$value=as.numeric(df$value)
    df$GTExTissue=as.factor(df$GTExTissue)
  
    # graph variant capture, prediction rate across all methods and tissues
    small.rates = filter(df, metric %in% c("% variants in \nany enhancer","Prediction rate"))
    small.rates$metric = plyr::revalue(small.rates$metric, 
                                       c("% variants in \nany enhancer"="Overlaps predicted \nenhancer",
                                         "Prediction rate"="Overlaps predicted enhancer \nlinked to correct gene"))
    sr = ggplot(small.rates, aes(x=metric, y=value, fill=method)) + 
      geom_bar(stat="identity", position="dodge", width=0.5) + 
      facet_grid(GTExTissue~., scales='fixed') + theme_minimal() + xlab('') + ylab('Fraction of fine-mapped, distal noncoding\neQTL variants with PIP >= 0.5') + 
      #theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=8)) +
      theme(text = element_text(size = rel(4)), axis.text.y=element_text(size=rel(2.5)),legend.text=element_text(size=rel(2.5))) +
      #scale_fill_viridis(discrete=TRUE,option='viridis',name='Method')
      scale_fill_manual(values=cpList,name='Method')
    
  
    pdf(file=outFile.pred,width=8,height=nrow(GTExTissues)*2+2); print(sr); dev.off()
    
    # graph positive predictive value
    df.ppv = filter(df, metric=="Positive predictive value")
    
    #reorder method factor for graphing
    df.ppv$method=as.factor(df.ppv$method); 
    ppv = ggplot(df.ppv, aes(x=method, y=value, fill=method)) + 
      geom_bar(stat="identity", position="dodge", width=0.5) + 
      facet_grid(GTExTissue~., scales='fixed') + theme_minimal() + xlab('') + ylab('Positive predictive value') + 
      theme(text = element_text(size = rel(4)),axis.text.x = element_text(angle=60,hjust=1), legend.position='none') + 
      scale_fill_manual(values=cpList)
    
    pdf(file=outFile.PPV,width=8,height=nrow(GTExTissues)*2+2); print(ppv); dev.off()
    
    
    # print table
    df$metric = str_remove(df$metric, "\n")
    write.table(df, file=paste0(outDir, "/predictionMetrics.tsv"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
    
}

main()
