main <- function() {
    suppressPackageStartupMessages({library(ggplot2)
                                    library(dplyr)
                                    library(scales)
                                    library(tidyr)
                                    library(tibble)
                                    library(viridis)
                                    library(optparse)
                                    library(egg)})
    
    # load data
    option_list <- list(
        make_option(c("--outfile"), type="character", default=NA, help="output filepath"),
        make_option(c("--useCategory"), type="character", default="False", help="if 'True' produce heatmap with category, if provided in biosample key; otherwise, don't"),
        make_option(c("--table"), type="character", default=NA, help="enrichment table"),
        make_option(c("--enhancersizes"), type="character", default=NA, help= "number of base pairs per enhancer set"),
        make_option(c("--samplekey"), type="character", default=NA, help=".csv file linking cell type IDs to sample names or category (optional)"),
        make_option(c("--maxcolor"), type="numeric", default=NA, help="max enrichment score to use for color scale")
    )
    
    opt <- parse_args(OptionParser(option_list=option_list))
    
    file.name = opt$table; enhancer.size = opt$enhancersizes;  
    sampleKeyFile = opt$samplekey; out.file = opt$outfile
    max.color=opt$maxcolor; useCat = opt$useCategory
    
    enrMatrix = read.table(file=file.name, header=TRUE, stringsAsFactors = FALSE) %>% 
        filter(enrichment!='NA')

    if (is.na(max.color) || max.color=='nan') {
        max.color = ceiling(quantile(enrMatrix$enrichment, 0.95))
    }
    
    ## add category/sample name
    if (is.na(sampleKeyFile)){
        enrMatrix$sampleCategory = enrMatrix$sampleName
    } else {
        cat.data = read.table(sampleKeyFile, header=TRUE, sep="\t", fill=TRUE)
        if ("sampleCategory" %in% colnames(cat.data) && useCat=="True"){
            cat.data = dplyr::select(cat.data, biosample, sampleCategory)
            enrMatrix = left_join(enrMatrix, cat.data, by=c('Biosample'='biosample'))
        } else {
            enrMatrix$sampleCategory = enrMatrix$sampleName
        }
    }
    # add base pairs per  biosample
    bp.data = read.table(enhancer.size, header=FALSE)
    colnames(bp.data) = c('Biosample','EnhancerMb')
    bp.data$EnhancerMb = bp.data$EnhancerMb/1E6
    enrMatrix = left_join(enrMatrix, bp.data)

    # calc aggregate statistics
    cat.key = aggregate(cbind(EnhancerMb, enrichment, nVariantsOverlappingEnhancers) ~ sampleCategory + GTExTissue, data=enrMatrix, FUN=mean)
    colnames(cat.key)[3:5] = c('EnhancerMbByCategory', 'enrichmentByCategory', 'nVariantsOverlappingEnhancersByCategory')
    enrMatrix = left_join(enrMatrix, cat.key, by=c('sampleCategory', 'GTExTissue'))

    # cluster to get orders
    M = dplyr::select(enrMatrix, sampleCategory, GTExTissue, enrichmentByCategory) %>% distinct() %>% spread(GTExTissue, enrichmentByCategory) %>% column_to_rownames("sampleCategory") %>% drop_na()
    ord.GTEx = hclust(dist(1-cor(M)), method = "ward.D")$order
    ord.biosample = hclust(dist(1-cor(t(M))), method = "ward.D")$order
    enrMatrix$sampleCategory = factor(enrMatrix$sampleCategory)
    enrMatrix$sampleCategory = factor(enrMatrix$sampleCategory, levels(enrMatrix$sampleCategory)[ord.biosample])
    enrMatrix$GTExTissue = factor(enrMatrix$GTExTissue)
    enrMatrix$GTExTissue = factor(enrMatrix$GTExTissue, levels(enrMatrix$GTExTissue)[ord.GTEx])

    # heat map
    g = ggplot(data=enrMatrix,aes(x=GTExTissue,y=sampleCategory,fill=enrichmentByCategory)) + geom_tile()  +  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.text=element_text(size=12), legend.position='top', legend.text=element_text(size=12), legend.title=element_text(size=12)) + ylab('') + scale_fill_viridis(discrete = FALSE, name='Enrichment', limits=c(1,max.color),oob=scales::squish)
    
    g.solo = ggplot(data=enrMatrix,aes(x=GTExTissue,y=sampleCategory,fill=enrichmentByCategory)) + geom_tile()  +  theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 60,hjust=1), axis.text=element_text(size=12), legend.position='top', legend.text=element_text(size=12), legend.title=element_text(size=12)) + ylab('') + scale_fill_viridis(discrete = FALSE, name='Enrichment', limits=c(1,max.color),oob=scales::squish)
    
    
    # variants per tissue
    tissueCount = dplyr::select(enrMatrix,GTExTissue,nVariantsGTExTissue) %>% distinct()
    t = ggplot(tissueCount, aes(x=GTExTissue, y=nVariantsGTExTissue)) + geom_bar(stat='identity', width=0.5) + theme_minimal() + ylab('# variants in \n GTEx tissue') + theme(axis.text.x = element_text(angle = 60,hjust=1), axis.text=element_text(size=12), axis.title=element_text(size=12))  + xlab('')
    
    # Mb per enhancer set, all  biosamples
    e = ggplot(enrMatrix, aes(x=sampleCategory, y=EnhancerMbByCategory)) + geom_bar(stat='identity', width=0.5) + theme_minimal() + theme(axis.text=element_text(size=12), axis.title=element_text(size=12),  axis.text.y=element_blank()) + ylab('Enhancer\nset size (Mb)') + xlab('') + coord_flip()
    
    blank = ggplot() + theme_void()
    
    y=18; if (nlevels(enrMatrix$sampleCategory)>100) {y=nlevels(enrMatrix$sampleCategory)/4.5} # set output file size
    
    p3 = egg::ggarrange(g, e, t, blank, nrow=2, ncol=2, heights=c(0.9, 0.1),widths=c(0.9,.1))
    #p3 = g.solo
    pdf(out.file, width=20, height=y); print(p3); dev.off()

}

main()
