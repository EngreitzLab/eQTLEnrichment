suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(optparse)
  library(stringr)
  library(ggplot2)
  library(viridis)
  library(plyr)})

main <- function() {
  option_list <- list(
    make_option(c("--tables"), type="character", default=NA, help="tsv with variants"),
    make_option(c("--outDir"), type="character", default=NA, help="file path to results folder"))
  opt = parse_args(OptionParser(option_list=option_list))
  outDir = opt$outDir
  tables = strsplit(opt$tables, " ") %>% unlist()
  
  # get list of GTEx tissues represented
  GTExTissues = str_split(tables, '/') %>% data.frame() %>% t()
  GTExTissues = GTExTissues[,ncol(GTExTissues)] %>% str_split('\\.') %>% data.frame() %>% t() %>% data.frame() %>% dplyr::select(X1) %>% distinct()
  # add rates for proximity-based E-G predictions
  #df = data.frame(matrix(ncol=4,nrow=0, dimnames=list(NULL)))
  df = data.frame('temp','temp','temp',1)
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
  df = filter(df, GTExTissue!='temp')
  print('Check1')

    # read in prediction rates from each method
  for (fileName in tables){
    predTable = read.table(file=fileName, header=TRUE, stringsAsFactors=FALSE)
    tissue = str_split(fileName, '/')[[1]]; tissue = tissue[length(tissue)] %>% str_split(pattern='\\.'); 
    tissue = tissue[[1]][1]
    method = str_split(fileName, '/')[[1]]; method = method[length(method)-1]
    prediction.rate.GivenEnhancer = nrow(filter(predTable, predictionClass=='inEnhancer-correctGene'))/nrow(filter(predTable,predictionClass!='noOverlap'))
    prediction.rate.total = nrow(filter(predTable, predictionClass=='inEnhancer-correctGene'))/nrow(predTable)
    prediction.rate.inEnhancer=nrow(filter(predTable,predictionClass!='noOverlap'))/nrow(predTable)
    print("Check0")
    print(df)
    df = add_row(df, GTExTissue=tissue, metric="Prediction rate given\nvariant in enhancer", method=method, value=prediction.rate.GivenEnhancer)
    print("Check0.1")
    print(df)
    df = add_row(df, GTExTissue=tissue, metric="% variants in\nany enhancer", method=method, value=prediction.rate.inEnhancer)
    print("Check0.2")
    df = add_row(df, GTExTissue=tissue, metric="Prediction rate", method=method, value=prediction.rate.total)    
    
    #df = rbind(df, c(tissue, "Prediction rate given\nvariant in enhancer", method, prediction.rate.GivenEnhancer)) %>%
      #rbind(c(tissue, "% variants in\nany enhancer", method, prediction.rate.inEnhancer)) %>%
      #rbind(c(tissue, "Prediction rate", method, prediction.rate.total))
  }
  
  print('Check2')
  
  
  df$value=as.numeric(df$value)
  df$GTExTissue=as.factor(df$GTExTissue)
  df$metric = factor(df$metric); df$metric=factor(df$metric, levels = levels(df$metric)[c(2,3,6,1,4,5)])
  
  #### graph variant prediction rates for each set
  methods = df$method %>% data.frame() %>% distinct() %>% dplyr::filter(.!="Proximity") %>% distinct()

  for (method.i in methods[[1]]){
    for (tissue.i in GTExTissues[[1]]){
      #df.specific = data.frame(matrix(ncol=4,nrow=0, dimnames=list(NULL)))
      df.specific = data.frame('temp','temp','temp',1)
      key = paste0(method.i, "/", tissue.i)
      #temp = tables[str_detect(tables, method.i)]; temp = temp[str_detect(temp, tissue.i)]
      temp = tables[str_detect(tables, key)] 
      if (nrow(as.data.frame(temp))>0) { # if this method has predictions for this tissue
        # get list of all tables for this tissue
        tables.tissue = tables[str_detect(tables, paste0(tissue.i))]
        key.table = temp
        # filter all to list of variants in enhancers for this method
        var.set = read.table(file=key.table, header=TRUE, stringsAsFactors=FALSE) %>% filter(predictionClass!="noOverlap") %>% dplyr::select(chr, start, end) %>% distinct()
        # add prediction rates for each method 
        for (x in tables.tissue){
          t = read.table(file=x, header=TRUE, stringsAsFactors=FALSE)

          this.method = str_split(x, '/')[[1]]; this.method = this.method[length(this.method)-1]
          t = suppressMessages(left_join(var.set, t))
          rate = nrow(filter(t, predictionClass=='inEnhancer-correctGene'))/nrow(t)
          metric.name = paste0("Prediction rate given variant\nin enhancer (", this.method, ")" )
          colnames(df.specific) = c("GTExTissue", "metric", "method","value")
          df.specific = add_row(df.specific, GTExTissue=tissue.i, metric=metric.name, method=this.method, value=rate)
          #df.specific = rbind(df.specific, c(tissue.i, metric.name, this.method, rate))
          print('Check3')
        }
        df.specific = filter(df.specific, GTExTissue!='temp')
    
        # add proximity rates for this tissue & this method
        t = read.table(file=key.table, header=TRUE, stringsAsFactors=FALSE)
        t = suppressMessages(left_join(var.set, t)); 
        
        closest.gene.rate = nrow(filter(t, closestGene==eGene))/nrow(t)
        closest.TSS.rate = nrow(filter(t, closestTSS==eGene))/nrow(t)
        nearby.TSS.rate = nrow(filter(t, nearbyTSS==TRUE))/nrow(t)
        
        df.specific = add_row(df.specific, GTExTissue=tissue.i, metric="Closest gene body", method="Proximity", value=closest.gene.rate) %>% add_row(GTExTissue=tissue.i, metric="Closest TSS", method="Proximity", value=closest.TSS.rate) %>% add_row(GTExTissue=tissue.i, metric="Any TSS within 100 kb", method="Proximity", value=nearby.TSS.rate)
        #df.specific = rbind(df.specific, c(tissue.i, "Closest gene body","Proximity", closest.gene.rate)) %>% rbind(c(tissue.i, "Closest TSS", "Proximity", closest.TSS.rate)) %>% rbind(c(tissue.i, "Any TSS within 100 kb", "Proximity", nearby.TSS.rate))
        print('Check4')
        
        # graph
        df.specific$value=as.numeric(df.specific$value)
        nVar = nrow(var.set); 
        graph.title = paste0('eGene prediction rates for ', tissue.i, ' variants\noverlapping ', method.i, ' enhancers (n=', nVar, ")")
        g.specific = ggplot(df.specific, aes(x=metric, y=value, fill=method)) + 
          geom_bar(stat="identity", position="dodge", width=0.5) + ylab('Sensitivity') +
          theme_minimal() +  ylim(0,1) + ggtitle(graph.title) + scale_fill_viridis(discrete = TRUE) +
          theme(axis.text.x = element_text(angle = 60,hjust=1), axis.title.x = element_blank(), 
                text = element_text(size=14),
                plot.title= element_text(size=14),
                legend.position='none')
        
        file.name = paste0(outDir, "/", tissue.i, ".", method.i, "Enhancers.variantOverlapSensitivity.pdf")
        pdf(file=file.name,width=8,height=5); print(g.specific); dev.off()
      
      }
    }
  }
}
main()