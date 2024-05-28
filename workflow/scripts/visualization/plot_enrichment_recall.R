## NEW RECALL
## load libraries
library(dplyr)
library(ggplot2)
library(colorspace)
library(data.table)
library(tidyr)

## INPUTS

files = snakemake@input$enrichmentRecall_files %>% strsplit(" ") %>% unlist()
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
df = dplyr::left_join(df, cp, by="method")
df$key = paste0(df$pred_name_long, " (", df$Biosample, ")")
df[is.na(df)] = 0
df = df[order(df$threshold),]

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

df_binary = dplyr::filter(df_plot, nPoints==2, threshold==1)
df_plot = dplyr::filter(df_plot, nPoints!=2)

# plotting params
ylim = 50
pred_colors = cp$hex
names(pred_colors) = cp$key

g=ggplot(data=df_plot, aes(x=recall.linking, y=enrichment, color=key)) +
  geom_line(linewidth=1) +
  geom_point(data=df_binary, aes(x=recall.linking, y=enrichment, color=key), size=4) +
  scale_color_manual(values=pred_colors) +
  ylab("Enrichment (eQTLs vs. common variants)") + xlab("Recall (variants overlapping prediction linked to eGene)") +
  labs(col="Predictor") +
  coord_cartesian(ylim=c(0,ylim)) +
  theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), legend.text = element_text(size=7), legend.title=element_text(size=8), legend.position="bottom", legend.direction="vertical")

ggsave(out_plot, g, width=5, height=5)
write.table(df, out_table, quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")

