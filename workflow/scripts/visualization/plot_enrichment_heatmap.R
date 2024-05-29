suppressPackageStartupMessages({library(ggplot2)
        library(dplyr)
        library(scales)
        library(tidyr)
    	library(tibble)
        library(egg)
		library(data.table)})
    

main <- function() {
	# input data
	enrTableFile = snakemake@input$enrichmentTable # 0-30000kb table
	enhSizeFile = snakemake@input$enhancerSizes # per method; biosample / base pairs
	mapFile= snakemake@input$map
	outFile_combined = snakemake@output$outFile_combined
	outFile_solo = snakemake@output$outFile_alone

	enr = fread(enrTableFile, sep="\t", header=TRUE) %>% drop_na()
	enr$Biosample[enr$Biosample=="Cells_EBV-transformed_lymphocytes"] = "Cells_EBV_transformed_lymphocytes"
	enhSizes = fread(enhSizeFile, sep="\t", header=TRUE) 
	colnames(enhSizes) = c("Biosample", "enhBp") 
	map = fread(mapFile, sep="\t", header=TRUE)

    # add base pairs per biosample
	enr = left_join(enr, enhSizes, by="Biosample")
	enr$enhMb = enr$enhBp/1e6
	
    # cluster to get orders
    M = dplyr::select(enr, Biosample, GTExTissue, enrichment) %>% distinct() %>%
		pivot_wider(names_from=GTExTissue, values_from = enrichment) %>% column_to_rownames("Biosample") %>% drop_na()
	order_tissues =  hclust(dist(1-cor(M)), method = "ward.D")$order
	order_biosamples = hclust(dist(1-cor(t(M))), method = "ward.D")$order
	enr$Biosample = factor(enr$Biosample, levels=rownames(M)[order_biosamples], ordered=TRUE)
	enr$GTExTissue = factor(enr$GTExTissue, levels=colnames(M)[order_tissues], ordered=TRUE)

	# set plotting params
    #colors = c("#c5373d", "#f7f7f7", "#006eae") # red-white-blue
	colors = c("#f6eff7","#bdc9e1", "#67a9cf","#1c9099", "#016c59")
	na_color = "#ffffff"
	max_value = round(quantile(enr$enrichment, 0.9), 1) # 90th percentile enrichment
	lims = c(0, max_value) # or 0,1? also try not log?
	ht = ifelse(length(rownames(M))>50, 16, 8)

	# mark intersections labeled as "matches"
	enr$key = paste0(enr$GTExTissue, ".", enr$Biosample)
	map$key = paste0(map$tissue, ".", map$biosample)
	enr = mutate(enr, label = ifelse(key %in% map$key, "*", ""))

    # heat map alone
	just_enr = ggplot(enr, aes(x=GTExTissue, y=Biosample, fill=enrichment)) + 
		geom_tile() +
		geom_text(aes(label = label), size=6) +
		scale_fill_gradientn(colors=colors, oob=scales::squish, na.value="#FFFFFF", limits=lims, name="Enrichment") +
		theme_minimal() + theme(axis.text = element_text(size = 7), axis.title = element_blank(), axis.text.x = element_text(angle=60, hjust=1),
			legend.position='top',  legend.direction='horizontal', legend.text=element_text(size=7), legend.title=element_text(size=7))
		
	# plots for grid
	enr_grid  = ggplot(enr, aes(x=GTExTissue, y=Biosample, fill=enrichment)) + 
		geom_tile() +
		geom_text(aes(label = label), size=6) +
		scale_fill_gradientn(colors=colors, oob=scales::squish, na.value="#FFFFFF", limits=lims, name="Enrichment") +
		theme_minimal() + theme(axis.text = element_text(size = 7), axis.title = element_blank(), axis.text.x = element_blank(),
			legend.position='top',  legend.direction='horizontal', legend.text=element_text(size=7), legend.title=element_text(size=7))

	nVar = dplyr::select(enr, GTExTissue, nVariantsGTExTissue) %>% distinct()
	var_count = ggplot(nVar, aes(x=GTExTissue, y=nVariantsGTExTissue)) +
		geom_bar(stat="identity", width=0.5) +
		ylab("# variants in tissue") +  xlab("") +
		theme_minimal() + theme(axis.text = element_text(size = 7), axis.text.x = element_text(angle=60, hjust=1))

	enhSizes = dplyr::select(enr, Biosample, enhMb) %>% distinct()
	enh_size = ggplot(enhSizes, aes(x=Biosample, y=enhMb)) +
		geom_bar(stat="identity", width=0.5) +
		ylab("Enhancer set size\n(Mb)") + xlab("") +
		theme_minimal() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), axis.text.y = element_blank()) + 
		coord_flip()

	blank = ggplot() + theme_void()

	assembled = egg::ggarrange(enr_grid, enh_size, var_count, blank, nrow=2, ncol=2, heights=c(2, 0.2), widths=c(2, 0.3))

	ggsave(outFile_solo, just_enr, width=8, height=ht)
	ggsave(outFile_combined, assembled, width=10, height=ht)

}

main()
