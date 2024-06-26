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
predTable_files = snakemake@input$predTable_files  %>% strsplit(" ") %>% unlist()
cp = fread(snakemake@input$colorPalette, sep="\t") # method, pred_name_long, hex
distances_min = snakemake@params$distances_min %>% as.character() %>% strsplit(" ") %>% unlist() %>% as.numeric()
distances_max = snakemake@params$distances_max %>% as.character() %>% strsplit(" ") %>% unlist() %>% as.numeric()
out_plot = snakemake@output$outFile


### FORMAT DATA
# aggregate enrichment tables
for (i in 1:length(enrTable_files)){
	enr = fread(file = enrTable_files[i], header = TRUE, sep="\t")  
	if(i==1 ){enr.all = enr
	} else {enr.all = rbind(enr.all, enr)}
}

# aggregate recall tables
for (i in 1:length(predTable_files)){
	rec = fread(file = predTable_files[i], header = TRUE, sep="\t")
	if(i==1) {rec.all = rec
	} else {rec.all = rbind(rec.all, rec)}
}

df = inner_join(enr.all, rec.all, by=c("distance_min", "distance_max", "method", "GTExTissue", "Biosample")) %>%
	dplyr::filter(!is.na(enrichment), !is.na(recall.linking)) %>%
	dplyr::filter(enrichment>0, recall.linking>0, total.variants>20)
df$log10_enrichment = log10(df$enrichment)

## make y-axis labels (distance range, min-max variants)
dist_labels = data.frame(distance_min=distances_min, distance_max=distances_max)
dist_labels$distance.label = ""
for (i in 1:nrow(dist_labels)){
	df.this = dplyr::filter(df, distance_max==dist_labels$distance_max[i])
	n_min = min(df.this$total.variants)
	n_max = max(df.this$total.variants)
	count = paste0("(N = ", n_min, "-", n_max, ")")
	if (dist_labels$distance_max[i]==30000) {
		cat = "All variants"
	} else {
		cat = paste0(dist_labels$distance_min[i],  "-", dist_labels$distance_max[i], " Kb")
	}
	dist_labels$count[i]  = count
	dist_labels$cat[i] = cat
	dist_labels$distance.label[i] = paste0(cat, "\n", count)
}
df = left_join(df, dist_labels, by=c("distance_min", "distance_max"))
df = left_join(df, cp, by="method")
df$pred_name_long = factor(df$pred_name_long, levels=cp$pred_name_long, ordered=TRUE)

## color palette
pred_colors = cp$hex
names(pred_colors) = cp$pred_name_long

### GENERATE PLOTS
plots =  vector(mode="list", length(distances_max))
enr_max = max(df$log10_enrichment)

## small multiples
g = ggplot(df, aes(x=recall.linking, y=log10_enrichment, color=pred_name_long)) +
	geom_point(alpha=0.5) +
	xlab("Recall (fraction of variants overlapping variant linked to eGene)") + ylab("log10 enrichment (eQTLs vs. common variants)") +
	scale_color_manual(values=pred_colors) +
	facet_grid(rows=vars(distance.label), cols=vars(pred_name_long)) +
	theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), aspect.ratio = 1, legend.position='None')

## just all dist
# df.all = dplyr::filter(df, distance_max==30000)
# x_label = paste0("Recall (fraction of variants overlapping variant linked to eGene)\n", df.all$count[1])
# g = ggplot(df.all, aes(x=recall.linking, y=log10_enrichment, color=pred_name_long)) +
# 	geom_point(alpha=0.5) +
# 	xlab(x_label) + ylab("log10 enrichment (eQTLs vs. common variants)") +
# 	scale_color_manual(values=pred_colors, name="Predictor") +
# 	theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), legend.position='right')
# h=5
# w = 8


h = nrow(dist_labels) * 1.5
w = nrow(cp) * 1.5
ggsave(out_plot, g, width=w, height=h)

