suppressPackageStartupMessages({
  library(plyr)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(egg)
  library(forcats)
  library(data.table)
  library(stringr)
})


### INPUTS
enrTable_files = snakemake@input$enrichmentTable_files  %>% strsplit(" ") %>% unlist()
cp = fread(snakemake@input$colorPalette, sep="\t") # method, pred_name_long, hex
distances_min = snakemake@params$distances_min %>% as.character() %>% strsplit(" ") %>% unlist() %>% as.numeric()
distances_max = snakemake@params$distances_max %>% as.character() %>% strsplit(" ") %>% unlist() %>% as.numeric()
out_plot = snakemake@output$outFile


### FORMAT DATA
## aggregate and process enrichment tables
for (i in 1:length(enrTable_files)){
	temp = fread(file = enrTable_files[i], header = TRUE, sep="\t")  
	temp$tissue.biosample  = paste0(temp$GTExTissue, ".", temp$Biosample)
	method.this = temp$method[1]

	if(i==1){enr.all = temp} else {
		enr.all = rbind(enr.all, temp)
	}
}
enr.all = dplyr::filter(enr.all, is.finite(enrichment))

## make y-axis labels (distance range, min-max variants)
dist_labels = data.frame(distance_min=distances_min, distance_max=distances_max)
dist_labels$distance.label = " "
for (i in 1:nrow(dist_labels)){
	enr.this = dplyr::filter(enr.all, distance_max==dist_labels$distance_max[i])
	n_min = min(enr.this$nVariantsGTExTissue)
	n_max = max(enr.this$nVariantsGTExTissue)
	count = paste0("N = ", n_min, "-", n_max, " variants")
	if (dist_labels$distance_max[i]==30000) {
		cat = "All variants"
	} else {
		cat = paste0(dist_labels$distance_min[i],  "-", dist_labels$distance_max[i], " Kb")
	}
	dist_labels$count[i]  = count
	dist_labels$cat[i] = cat
	dist_labels$distance.label[i] = paste0(cat, "\n", count)
}
enr.all = left_join(enr.all, dist_labels, by=c("distance_min", "distance_max"))
enr.all = left_join(enr.all, cp, by="method")
enr.all$pred_name_long = factor(enr.all$pred_name_long, levels=cp$pred_name_long, ordered=TRUE)

## color palette
pred_colors = cp$hex
names(pred_colors) = cp$pred_name_long

# for biosamples_predictors plot 
# ordered.methods = c('In element (DHS) & closest gene', 'ABC_A=DNase, C=Avg. Intact Hi-C', 'ENCODE-E2G', 'EpiMap', 'EPIraction', 'ABC_A=DNase x H3K27ac, C=Avg. Intact Hi-C')
# enr.all$pred_name_long = factor(enr.all$pred_name_long, levels=ordered.methods)

### GENERATE PLOTS

g = ggplot(enr.all, aes(x=distance.label, y=enrichment, fill=pred_name_long)) +
	geom_violin(color="transparent", adjust = .5, scale="width") +
	#geom_point(position = position_jitterdodge(seed = 1, dodge.width = 0.7), size=0.5) +
	xlab("eVariant - eGene distance") + ylab("Enrichment (eQTLs vs. common variants)") +
	scale_fill_manual(values=pred_colors, name="Predictor") +
	theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8),
		legend.position='right', legend.direction='vertical', legend.title=element_text(size=8), legend.text=element_text(size=7))

ggsave(out_plot, g, width=8, height=4)

