## NEW RECALL
## load libraries
library(dplyr)
library(ggplot2)
library(colorspace)
library(data.table)
library(tidyr)

## INPUTS
files = snakemake@input$variantsPerTissue %>% strsplit(" ") %>% unlist()
methods = snakemake@params$methods %>% strsplit(" ") %>% unlist()
tissues_matched = snakemake@params$tissues_matched %>% strsplit(" ") %>% unlist()
out_plot = snakemake@output$out_plot
out_plot_matched = snakemake@output$out_plot_matched
out_table = snakemake@output$out_table
distances_min = snakemake@params$distances_min %>%  as.character() %>% strsplit(" ") %>% unlist() %>% as.numeric()
distances_max = snakemake@params$distances_max %>%  as.character() %>% strsplit(" ") %>% unlist() %>% as.numeric()

# gather and average tables
n = length(methods)

for (i in 1:n){
	temp = fread(files[i], sep="\t")
	temp = pivot_longer(temp, -tissue, names_to="distance_bin", values_to="var_count_temp")
	if (i==1){
		df = temp %>% mutate(var_count = var_count_temp) %>% dplyr::select(tissue, var_count, distance_bin)
		} else {
		df = left_join(df, temp, by=c("tissue", "distance_bin")) %>%
			mutate(var_count = var_count + var_count_temp) %>%
			dplyr::select(tissue, var_count, distance_bin)
	}
}
df$var_count = round(df$var_count/n)

dist = data.frame(distance_min = distances_min, distance_max = distances_max)
dist$distance_bin = ""
dist$distance_label = ""
for (i in 1:nrow(dist)){
	if (dist$distance_max[i] == 30000){
		dist$distance_bin[i] = "n_all"
		dist$distance_label[i] = "All variants"
	} else {
		dist$distance_bin[i] = paste0("n_bin", i)
		dist$distance_label[i] = paste0(dist$distance_min[i], "-", dist$distance_max[i], " Kb")
	}
}
last_row = c(distance_min=distances_max[length(distances_max)-1], distance_max=distances_max[length(distances_max)], distance_bin=paste0("n_bin", length(distances_max)), distance_label=paste0("> ", distances_max[length(distances_max)-1], "Kb"))
dist = rbind(dist, last_row)

dist = dplyr::filter(dist, distance_bin!="n_all")
colors_all = c("#003648", "#006479", "#0096a0", "#49bcbc", "#96ced3", "#cae5ee") # dark to light, assume you have <=6 distance bins
dist$hex = colors_all[1:nrow(dist)]
cp = dist$hex
names(cp) = dist$distance_label

to_order = dplyr::filter(df, distance_bin=="n_all") %>% arrange(var_count)

df = dplyr::filter(df, distance_bin!="n_all") %>%
	left_join(dist, by="distance_bin")
df$tissue = factor(df$tissue, levels=to_order$tissue, ordered=TRUE)
df$distance_label = factor(df$distance_label, levels=dist$distance_label, ordered=TRUE)

g=ggplot(data=df, aes(x=tissue, y=var_count, fill=distance_label)) +
	geom_hline(yintercept=50, linetype="dashed", color="#96a0b3") +
	geom_bar(position="stack", stat="identity") +
	labs(x="", y="Number of eQTLs above PIP threshold", fill="eVariant-eGene distance") +
	scale_fill_manual(values=cp) +
	theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8)) +
	coord_flip() 

n_tissues = length(unique(df$tissue))
ht = max(4, n_tissues/12)
ggsave(out_plot, g, width=6, height=ht)
write.table(df, out_table, quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")

# plot just matched tissues
df_filt = dplyr::filter(df, tissue %in% tissues_matched)
g=ggplot(data=df_filt, aes(x=tissue, y=var_count, fill=distance_label)) +
	#geom_hline(yintercept=50, linetype="dashed", color="#96a0b3") +
	geom_bar(position="stack", stat="identity") +
	labs(x="", y="Number of eQTLs above PIP threshold", fill="eVariant-eGene distance") +
	scale_fill_manual(values=cp) +
	theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8)) +
	coord_flip() 

n_tissues = length(unique(df_filt$tissue))
ht = max(4, n_tissues/12)
ggsave(out_plot_matched, g, width=6, height=ht)

