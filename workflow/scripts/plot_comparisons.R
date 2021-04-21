# take enrichment table produced by enrichment_matrix.sh, bp per enhancer, and output directory as input arg
# also adjust scale of heat map according to range of enrichment values

main <- function() {
  suppressPackageStartupMessages({library(ggplot2)
                                  library(dplyr)
                                  library(tidyr)
                                  library(optparse)
                                  library(viridis)})

    # load data
  
    enrichmentTables = (snakemake@input$enrichmentTables) %>% strsplit(" ") %>% unlist()
    names = (snakemake@params$methods) %>% strsplit(" ") %>% unlist()
    outCDF = (snakemake@output$cdf)
    outDensity = (snakemake@output$density)
    outBoxplot = (snakemake@output$boxplot)
  
    # aggregate data
    enr.all = read.table(file=enrichmentTables[1], header=TRUE, stringsAsFactors = FALSE)
    enr.all$predictionSet = names[1]
    
    if(length(names)>1){
      for (i in 2:length(names)){
        temp = read.csv(file=enrichmentTables[i], sep='\t', header=TRUE, stringsAsFactors = FALSE)
        temp$predictionSet = names[i]
        enr.all = rbind(enr.all, temp)
      }
    }
    
    enrLabel = 'Enrichment\n(GTEx variants/all common variants)'
    cdf = ggplot(enr.all, aes(enrichment, col=predictionSet)) + stat_ecdf(geom = "step") + 
      ylab('Cumulative fraction') + xlab(enrLabel) + 
      theme_minimal() + scale_color_viridis(discrete=TRUE,name='Method') + 
      theme(text = element_text(size = rel(4)), legend.text=element_text(size=rel(2.5)))
    
    d = ggplot(enr.all,aes(x=enrichment, col=predictionSet)) + geom_density() + 
      theme_minimal() + ylab('Density') + xlab(enrLabel) + 
      scale_color_viridis(discrete=TRUE,name='Method') + 
      theme(text = element_text(size = rel(4)), legend.text=element_text(size=rel(2.5)))
    
    bp = ggplot(enr.all, aes(x=predictionSet, y=enrichment, fill=predictionSet)) + geom_boxplot() +
      theme_minimal() + ylab(enrLabel) + xlab('') + 
      scale_fill_viridis(discrete=TRUE) +
      theme(text = element_text(size = rel(4)), legend.position='none', 
            axis.text.x=element_text(angle=60,hjust=1))
      
      
    pdf(file=outCDF, width=7, height=5); print(cdf); dev.off()
    pdf(file=outDensity, width=7, height=5); print(d); dev.off()
    pdf(file=outBoxplot, width=7, height=5); print(bp); dev.off()
    
    
    }

main()
