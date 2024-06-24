suppressPackageStartupMessages({library(ggplot2)
        library(dplyr)
        library(scales)
        library(tidyr)
    	library(tibble)
        library(egg)
		library(data.table)})
    

main <- function() {
	# input data
	recallTableFile = snakemake@input$recallTable # by distance
	enhSizeFile = snakemake@input$enhancerSizes # per method; biosample / base pairs
	p_threshold = snakemake@params$p_threshold %>% as.numeric()
	outFile_combined = snakemake@output$outFile_combined
	outFile_solo = snakemake@output$outFile_alone

	recall = fread(recallTableFile, sep="\t", header=TRUE)
	recall = dplyr::filter(recall, distance_max==30000, is.finite(recall.linking), !is.na(recall.linking), total.variants>20)
	# remove tissues where sd recall = 0
	recall_sum = group_by(recall, GTExTissue) %>%
		summarize(sd_recall = sd(recall.linking)) %>%
		dplyr::filter(sd_recall>0)
	recall = dplyr::filter(recall, GTExTissue %in% recall_sum$GTExTissue)

	recall$Biosample[recall$Biosample=="Cells_EBV-transformed_lymphocytes"] = "Cells_EBV_transformed_lymphocytes"
	enhSizes = fread(enhSizeFile, sep="\t", header=TRUE) 
	colnames(enhSizes) = c("Biosample", "enhBp") 

    # add base pairs per biosample
	recall = left_join(recall, enhSizes, by="Biosample")
	recall$enhMb = recall$enhBp/1e6
	
    # cluster to get orders
    M = dplyr::select(recall, Biosample, GTExTissue, recall.linking) %>% distinct() %>%
		pivot_wider(names_from=GTExTissue, values_from = recall.linking) %>% column_to_rownames("Biosample") %>% drop_na()
	print(M)
	print(cor(M))
	order_tissues =  hclust(dist(1-cor(M)), method = "ward.D2")$order
	order_biosamples = hclust(dist(1-cor(t(M))), method = "ward.D2")$order
	recall$Biosample = factor(recall$Biosample, levels=rownames(M)[order_biosamples], ordered=TRUE)
	recall$GTExTissue = factor(recall$GTExTissue, levels=colnames(M)[order_tissues], ordered=TRUE)

	# set plotting params
	colors = c("#edf8fb","#b3cde3", "#8c96c6", "#8856a7", "#810f7c")
	na_color = "#ffffff"
	lims = c(0, 0.25) 
	ht = ifelse(length(rownames(M))>50, 16, 8)

    # heat map alone
	just_recall = ggplot(recall, aes(x=GTExTissue, y=Biosample, fill=recall.linking)) + 
		geom_tile() +
		scale_fill_gradientn(colors=colors, oob=scales::squish, na.value="#FFFFFF", limits=lims, name="Recall (linked to correct eGene)") +
		theme_minimal() + theme(axis.text = element_text(size = 7), axis.title = element_blank(), axis.text.x = element_text(angle=60, hjust=1),
			legend.position='top',  legend.direction='horizontal', legend.text=element_text(size=7), legend.title=element_text(size=7))
		
	# plots for grid
	recall_grid  = ggplot(recall, aes(x=GTExTissue, y=Biosample, fill=recall.linking)) + 
		geom_tile() +
		scale_fill_gradientn(colors=colors, oob=scales::squish, na.value="#FFFFFF", limits=lims, name="Recall (linked to correct eGene)") +
		theme_minimal() + theme(axis.text = element_text(size = 7), axis.title = element_blank(), axis.text.x = element_blank(),
			legend.position='top',  legend.direction='horizontal', legend.text=element_text(size=7), legend.title=element_text(size=7))

	nVar = dplyr::select(recall, GTExTissue, total.variants) %>% distinct()
	var_count = ggplot(nVar, aes(x=GTExTissue, y=total.variants)) +
		geom_bar(stat="identity", width=0.5) +
		ylab("# eGene/eVariant pairs") +  xlab("") +
		theme_minimal() + theme(axis.text = element_text(size = 7), axis.text.x = element_text(angle=60, hjust=1))

	enhSizes = dplyr::select(recall, Biosample, enhMb) %>% distinct()
	enh_size = ggplot(enhSizes, aes(x=Biosample, y=enhMb)) +
		geom_bar(stat="identity", width=0.5) +
		ylab("Enhancer set size\n(Mb)") + xlab("") +
		theme_minimal() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), axis.text.y = element_blank()) + 
		coord_flip()

	blank = ggplot() + theme_void()

	assembled = egg::ggarrange(recall_grid, enh_size, var_count, blank, nrow=2, ncol=2, heights=c(2, 0.2), widths=c(2, 0.3))

	ggsave(outFile_solo, just_recall, width=8, height=ht)
	ggsave(outFile_combined, assembled, width=10, height=ht)

}

main()
