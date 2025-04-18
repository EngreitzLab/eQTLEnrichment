## NEW RECALL
## load libraries
library(dplyr)
library(ggplot2)
library(colorspace)
library(data.table)
library(tidyr)

## INPUTS

files = snakemake@input$enrichmentRecall_files %>% strsplit(" ") %>% unlist()
score_thresholds = snakemake@params$score_thresholds %>% as.character() %>% strsplit(" ") %>% unlist() %>% as.numeric()
method_names = snakemake@params$methods %>% strsplit(" ") %>% unlist()
this_tissue = snakemake@wildcards$GTExTissue
out_plot = snakemake@output$er_combined
out_table = snakemake@output$er_combined_table
cpFilePlotting = snakemake@input$colorPalette
cp = fread(cpFilePlotting, sep="\t") # method, pred_name_long, hex

# gather tables
for (i in 1:length(files)){
	temp = fread(files[i], sep="\t")
	temp$nPoints = nrow(temp)

	if (i==1){df = temp} else {
		df = rbind(df, temp)
	}
}

# organize for plotting
scores = data.frame(method=method_names, score_threshold=score_thresholds)
cp = left_join(cp, scores, by="method")

df = dplyr::left_join(df, cp, by="method")
df$key = paste0(df$pred_name_long, " (", df$Biosample, ")")
df[is.na(df)] = 0
df = df[order(df$threshold),]

df <- dplyr::filter(df, total.variants > 20, recall.total > 0.001)

# handle multiple matches
cp = dplyr::select(df, method, key, hex) %>% distinct() 
methods = unique(cp$method)
if (length(methods) < length(unique(cp$key))) {
	# set new hex codes
	for (i in 1:length(methods)){
		cp.this =  dplyr::filter(cp, method==methods[i])
		n = nrow(cp.this)
		if (n>1){
			hex_base = cp.this$hex[1]
			cols = colorRampPalette(c(hex_base, '#000000'))(n+2) # range of colors from hex_base to black (where black is 2 above the base)
			for (k in 1:nrow(cp.this)){
				cp$hex[cp$key==cp.this$key[k]] = cols[k]  # set corresponding hex code
			}
		}
	}
}

df = dplyr::select(df, -hex) %>% left_join(cp, by=c("key", "method")) # for saving later!

# data for binary predictors to be plotted as points
df$key = factor(df$key)
df_plot = dplyr::filter(df, recall.linking>0) # remove (0,0) points

df_binary = dplyr::filter(df_plot, nPoints==2) %>% dplyr::filter(threshold==1)
df_thresh = dplyr::filter(df, threshold==score_threshold) # suggested thresholds
df_plot = dplyr::filter(df_plot, nPoints>2)

# plotting params
ylim = min(50, max(df_plot$enrichment))
pred_colors = cp$hex
names(pred_colors) = cp$key
n_keys = cp$key %>% unique() %>% length(); print(n_keys)
n_legend_cols = ceiling(n_keys / 12); print(n_legend_cols)

if (this_tissue == "AllMatches") {
	n_var_label = paste0(mean(df_plot$total.variants), " variant-biosample pairs")
} else {
	n_var_label = paste0(mean(df_plot$total.variants), " variants")
}
x_label = paste0("Recall (variants overlapping prediction linked to eGene)\n", n_var_label)


g=ggplot(data=df_plot, aes(x=recall.linking, y=enrichment, color=key)) +
  geom_line(linewidth=0.75) +
  geom_point(data=df_binary, aes(x=recall.linking, y=enrichment, color=key), size=3, shape = 16) +
  geom_linerange(data=df_thresh, aes(ymin=CI_enr_low, ymax=CI_enr_high), linewidth = 0.75) +
  geom_point(data=df_thresh, aes(x=recall.linking, y=enrichment, color=key), size=3, shape = 16) +
  scale_color_manual(values=pred_colors) +
  ylab("Enrichment (eQTLs versus common variants)") + xlab(x_label) +
  labs(col="Predictor") + 
  coord_cartesian(ylim=c(0,ylim)) +
  theme_classic() + theme(axis.text = element_text(size = 7, color = "#000000"),
  	axis.title = element_text(size = 8), axis.ticks = element_line(color = "#000000"),
	legend.text = element_text(size= 7), legend.title=element_text(size=8), legend.position="right", legend.direction="vertical",
	aspect.ratio = 1) + 
  guides(col = guide_legend(nrow = 12))

width = 3 + 3 * n_legend_cols
ggsave(out_plot, g, width=width, height=4)
write.table(df, out_table, quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")

