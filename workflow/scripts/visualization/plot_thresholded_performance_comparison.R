suppressPackageStartupMessages({
  library(plyr)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(egg)
  library(colorspace)
  library(forcats)
  library(data.table)
  library(stringr)
})


### INPUTS
enrTable_files = snakemake@input$enrichmentTable_files  %>% strsplit(" ") %>% unlist()
predTable_files = snakemake@input$predTable_files  %>% strsplit(" ") %>% unlist()
map_files = snakemake@input$map_files  %>% strsplit(" ") %>% unlist()
cp = fread(snakemake@input$colorPalette, sep="\t") # method, pred_name_long, hex
distances_min = snakemake@params$distances_min %>% as.character() %>% strsplit(" ") %>% unlist() %>% as.numeric()
distances_max = snakemake@params$distances_max %>% as.character() %>% strsplit(" ") %>% unlist() %>% as.numeric()
out_plot = snakemake@output$outFile
out_enrTable= snakemake@output$enrAllTable
out_predMetrics = snakemake@output$predictionMetrics

### FORMAT DATA
## aggregate and process enrichment tables
maps = data.frame(file_name = map_files)
maps$method = ""
for (i in 1:nrow(maps)){
	temp = strsplit(maps$file_name[i], "/") %>% unlist()
	method = temp[length(temp)-2] # -1 = "int"
	maps$method[i] = method
}

for (i in 1:length(enrTable_files)){
	temp = fread(file = enrTable_files[i], header = TRUE, sep="\t")  
	temp$tissue.biosample  = paste0(temp$GTExTissue, ".", temp$Biosample)
	method.this = temp$method[1]
	map.this = fread(maps$file_name[maps$method==method.this], sep="\t")
	map.this$key = paste0(map.this$tissue, ".", map.this$biosample)
	temp = dplyr::filter(temp, tissue.biosample %in% map.this$key)

	if(i==1){enr.all = temp} else {
		enr.all = rbind(enr.all, temp)
	}
}

# get stats, edit names, define labels 
enr.all = enr.all[order(enr.all$enrichment), ]
max.enr = max(enr.all$enrichment[is.finite(enr.all$enrichment)])
enrLabel = 'Enrichment\n(GTEx variants/all common variants)'
enr.all = left_join(enr.all, cp, by="method")
enr.all$pred_name_long = factor(enr.all$pred_name_long, levels=cp$pred_name_long, ordered=TRUE)

## aggregate and process prediction metrics
for (i in 1:length(predTable_files)){
	temp = fread(file = predTable_files[i], header = TRUE, sep="\t")  
	if(i==1){pred.all = temp} else {
		pred.all = rbind(pred.all, temp)
	}
}
pred.all = left_join(pred.all, cp, by="method")
pred.all$pred_name_long = factor(pred.all$pred_name_long, levels=cp$pred_name_long, ordered=TRUE)

## make y-axis labels (distance range, min-max variants)
dist_labels = data.frame(distance_min=distances_min, distance_max=distances_max)
dist_labels$distance.label = " "
for (i in 1:nrow(dist_labels)){
	pred.this = dplyr::filter(pred.all, distance_max==dist_labels$distance_max[i])
	n_min = min(pred.this$total.variants)
	n_max = max(pred.this$total.variants)
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
pred.all = left_join(pred.all, dist_labels, by=c("distance_min", "distance_max"))
enr.all = left_join(enr.all, dist_labels, by=c("distance_min", "distance_max"))
pred_colors = cp$hex
names(pred_colors) = cp$pred_name_long

# for biosamples_predictors plot 
# ordered.methods = c('In element (DHS) & closest gene', 'ABC_A=DNase, C=Avg. Intact Hi-C', 'ENCODE-E2G', 'EpiMap', 'EPIraction', 'ABC_A=DNase x H3K27ac, C=Avg. Intact Hi-C')
# enr.all$pred_name_long = factor(enr.all$pred_name_long, levels=ordered.methods)
# pred.all$pred_name_long = factor(pred.all$pred_name_long, levels=ordered.methods)

### GENERATE PLOTS
## enrichment
enr.boxplot = ggplot(enr.all, aes(x = distance.label, y = enrichment, fill = pred_name_long)) +
  geom_boxplot(linewidth = 0.5) +
  coord_flip() +
  theme_minimal() + ylab("Enrichment\n(GTEx variants/all common variants)") + xlab('') +
  scale_fill_manual(values=pred_colors) +
  theme(legend.position = 'none') +
  ggtitle("Enrichment of variants\nin predicted enhancers")

## overlaps predicted enhancer : recall.total
sr.overlaps = ggplot(pred.all, aes(x = distance.label, y = recall.total, fill=pred_name_long)) +
  geom_boxplot(linewidth = 0.5) +
  scale_fill_manual(values=pred_colors) +
  theme_minimal() + ggtitle('Variants overlapping\npredicted enhancers') + 
  ylab('Fraction of GTEx variants') + xlab('') +
  theme(axis.text.y = element_blank(), legend.position='none') + coord_flip()

## linked to correct eGene
sr.predicted = ggplot(pred.all, aes(x = distance.label, y = correctGene.ifOverlap, fill=pred_name_long)) +
  geom_boxplot(linewidth = 0.5) +
  scale_fill_manual(values=pred_colors, name="Predictor") +
  theme_minimal() + ggtitle('Variants linked to correct gene,\ngiven overlapping predicted enhancer') + xlab('') +
  ylab('Fraction of GTEx variants\noverlapping predicted enhancers') +
  theme(axis.text.y = element_blank()) + coord_flip()

## save final plots
pdf(file=out_plot, width=12, height=5)
	all.tissues =  ggarrange(enr.boxplot, sr.overlaps, sr.predicted, nrow=1, ncol=3)
dev.off()

pred.all = dplyr::select(pred.all, -distance.label)
enr.all = dplyr::select(enr.all, -distance.label)

write.table(pred.all, out_predMetrics, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
write.table(enr.all, out_enrTable, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
