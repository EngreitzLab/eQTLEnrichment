# take enrichment table produced by enrichment_matrix.sh, bp per enhancer, and output directory as input arg
# also adjust scale of heat map according to range of enrichment values

main <- function() {
    library(ggplot2); library(dplyr); library(scales); library(tidyr); library(tibble); library(viridis);
    library(optparse)
    
    # load data
    option_list <- list(
        make_option(c("--outfile"), type="character", default=NA, help="output filepath"),
        make_option(c("--table"), type="character", default=NA, help="enrichment table"),
        make_option(c("--enhancersizes"), type="character", default=NA, help= "number of base pairs per enhancer set"),
        make_option(c("--samplekey"), type="character", default=NA, help=".csv file linking cell type IDs to sample names or category (optional)"),
        make_option(c("--cellid"), type="character", default=NA, help="name of column with cell type ID (optional)"),
        make_option(c("--category"), type="character", default=NA, help="name of column with relevant category or sample name (optional)"),
        make_option(c("--maxcolor"), type="numeric", default=NA, help="max enrichment score to use for color scale")
    )
    
    opt <- parse_args(OptionParser(option_list=option_list))
    
    file.name = opt$table; enhancer.size = opt$enhancersizes;  sample.key = opt$samplekey; out.file = opt$outfile
    cellID = opt$cellid; sample.name = opt$category; max.color=opt$maxcolor
    
    enrMatrix = read.csv(file=file.name, sep='\t', header=TRUE, stringsAsFactors = FALSE) %>% filter(enrichment!='NA')

    if (is.na(max.color) || max.color=='nan') {
        max.color = ceiling(quantile(enrMatrix$enrichment, 0.9))
    }
    
    # add category/sample name
    if (is.na(sample.name) || sample.name=="NA" || sample.name=="nan"){
        enrMatrix$ABCCategory = enrMatrix$ABCBiosample
    } else {
        ABC.cat = read.csv(sample.key, sep=',', header=TRUE) %>% dplyr::select(cellID, sample.name);
        colnames(ABC.cat) = c('ABCBiosample','ABCCategory') 
        enrMatrix = left_join(enrMatrix, ABC.cat) %>% drop_na()
    }
    
    # add base pairs per ABC biosample
    ABC.data = read.csv(enhancer.size, sep='\t', header=FALSE)
    colnames(ABC.data) = c('ABCBiosample','ABCEnhancerMb')
    ABC.data$ABCEnhancerMb = ABC.data$ABCEnhancerMb/1E6
    enrMatrix = left_join(enrMatrix, ABC.data) %>% filter(!is.na(ABCCategory))
    
    # calc aggregate statistics
    cat.key = aggregate(cbind(ABCEnhancerMb, enrichment, nVariantsOverlappingABCEnhancers) ~ ABCCategory + GTExTissue, data=enrMatrix, FUN=mean)
    colnames(cat.key)[3:5] = c('ABCEnhancerMbByCategory', 'enrichmentByCategory', 'nVariantsOverlappingABCEnhancersByCategory')
    enrMatrix = left_join(enrMatrix, cat.key, by=c('ABCCategory', 'GTExTissue'))
    
    print('check 3')

    # cluster to get orders
    M = dplyr::select(enrMatrix, ABCCategory, GTExTissue, enrichmentByCategory) %>% distinct() %>% spread(GTExTissue, enrichmentByCategory) %>% column_to_rownames("ABCCategory") %>% drop_na()
    ord.GTEx = hclust(dist(1-cor(M)), method = "ward.D")$order
    ord.ABC = hclust(dist(1-cor(t(M))), method = "ward.D")$order
    enrMatrix$ABCCategory = factor(enrMatrix$ABCCategory)
    enrMatrix$ABCCategory = factor(enrMatrix$ABCCategory, levels(enrMatrix$ABCCategory)[ord.ABC])
    #levels(enrMatrix$ABCCategory)
    enrMatrix$GTExTissue = factor(enrMatrix$GTExTissue)
    enrMatrix$GTExTissue = factor(enrMatrix$GTExTissue, levels(enrMatrix$GTExTissue)[ord.GTEx])
    #levels(enrMatrix$GTExTissue)
    
    # heat map
    g = ggplot(data=enrMatrix,aes(x=GTExTissue,y=ABCCategory,fill=enrichmentByCategory)) + geom_tile()  +  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.text=element_text(size=12), legend.position='top', legend.text=element_text(size=12), legend.title=element_text(size=12)) + ylab('') + scale_fill_viridis(discrete = FALSE, name='Enrichment', limits=c(1,max.color),oob=scales::squish)
    
    g.solo = ggplot(data=enrMatrix,aes(x=GTExTissue,y=ABCCategory,fill=enrichmentByCategory)) + geom_tile()  +  theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 60,hjust=1), axis.text=element_text(size=12), legend.position='top', legend.text=element_text(size=12), legend.title=element_text(size=12)) + ylab('') + scale_fill_viridis(discrete = FALSE, name='Enrichment', limits=c(1,max.color),oob=scales::squish)
    
    
    # variants per tissue
    tissueCount = dplyr::select(enrMatrix,GTExTissue,nVariantsGTExTissue) %>% distinct()
    t = ggplot(tissueCount, aes(x=GTExTissue, y=nVariantsGTExTissue)) + geom_bar(stat='identity', width=0.5) + theme_minimal() + ylab('# variants in \n GTEx tissue') + theme(axis.text.x = element_text(angle = 60,hjust=1), axis.text=element_text(size=12), axis.title=element_text(size=12))  + xlab('')
    
    # Mb per enhancer set, all ABC biosamples
    e = ggplot(ABC.data, aes(x=ABCCategory, y=ABCEnhancerMbByCategory)) + geom_bar(stat='identity', width=0.5) + theme_minimal() + theme(axis.text=element_text(size=12), axis.title=element_text(size=12),  axis.text.y=element_blank()) + ylab('Enhancer\nset size (Mb)') + xlab('') + coord_flip()
    
    blank = ggplot() + theme_void()
    
    y=18; if (nlevels(enrMatrix$ABCCategory)>140) {y=28} # set output file size
    
    pdf(out.file, width=20, height=y); g.solo; dev.off()
    #p3 = egg::ggarrange(g, e, t, blank, nrow=2, ncol=2, heights=c(0.9, 0.1),widths=c(0.9,.1))

}

main()
