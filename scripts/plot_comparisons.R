# take enrichment table produced by enrichment_matrix.sh, bp per enhancer, and output directory as input arg
# also adjust scale of heat map according to range of enrichment values

main <- function() {
    library(ggplot2); library(dplyr); library(tidyr);library(optparse)

    # load data
    option_list <- list(
        make_option(c("--names"), type="character", default=NA, help="list of prediction names"),
        make_option(c("--tables"), type="character", default=NA, help="list of enrichment tables"),
        make_option(c("--outdir"), type="character", default=NA, help="out directory")
    )
  

    opt <- parse_args(OptionParser(option_list=option_list))
    #names = read.csv(opt$names, header=FALSE, sep=' ')
    #enrichmentTables = read.csv(opt$tables, header=FALSE, sep=' ')
    names = strsplit(opt$names, " ") %>% unlist()
    print(names)
    enrichmentTables = strsplit(opt$tables, " ") %>% unlist()
    outDir = opt$outdir
    
    
    
    # aggregate data
    enr.all = read.csv(file=enrichmentTables[1], sep='\t', header=TRUE, stringsAsFactors = FALSE)
    enr.all$predictionSet = names[1]
    
    for (i in 2:length(names)){
        temp = read.csv(file=enrichmentTables[i], sep='\t', header=TRUE, stringsAsFactors = FALSE)
        temp$predictionSet = names[i]
        enr.all = rbind(enr.all, temp)
    }
    
    cdf = ggplot(enr.all, aes(enrichment, col=predictionSet)) + stat_ecdf(geom = "step") + ylab('Cumulative fraction') + xlab('Enrichment (GTEx variants/all common variants)') + theme_minimal() + scale_color_discrete(name='Prediction set') + theme(text = element_text(size = rel(4)), legend.text=element_text(size=rel(2.5)))
    d = ggplot(enr.all,aes(x=enrichment, col=predictionSet)) + geom_density() + theme_minimal() + ylab('Density') + xlab('Enrichment (GTEx variants/all common variants)') + scale_color_discrete(name='Prediction set') + theme(text = element_text(size = rel(4)), legend.text=element_text(size=rel(2.5)))
    
    pdf(file=paste0(outDir, "/cdf.pdf"), width=7, height=5); print(cdf); dev.off()
    pdf(file=paste0(outDir, "/density.pdf"), width=7, height=5); print(d); dev.off()
    
    }

main()
