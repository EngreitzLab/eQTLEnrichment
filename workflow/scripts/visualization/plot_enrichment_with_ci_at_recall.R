## NEW RECALL
## load libraries
library(dplyr)
library(ggplot2)
library(colorspace)
library(data.table)
library(tidyr)
library(ggpubr)

## INPUTS
files = snakemake@input$enrichmentRecall_files %>% strsplit(" ") %>% unlist()
out_plot = snakemake@output$enr_at_recall
out_table = snakemake@output$enr_at_recall_table
out_sign = snakemake@output$sign_table
cpFilePlotting = snakemake@input$colorPalette
cp = fread(cpFilePlotting, sep="\t") # method, pred_name_long, hex
recall.this = snakemake@wildcards$recall %>% as.numeric()
sign_threshold = snakemake@params$thresholdPval %>% as.numeric()

# gather tables
for (i in 1:length(files)){
	temp = fread(files[i], sep="\t")
	temp = drop_na(temp)
	temp = dplyr::filter(temp,is.finite(recall.linking))

	if (i==1){
		df = data.frame(matrix(nrow = 0, ncol = length(colnames(temp)))) 
		colnames(df) = colnames(temp)
	}
	
	# is closest recal within 0.02 of input?
	n_thresholds = nrow(temp)
  	recall_low = min(temp$recall.linking)
  	recall_high = max(temp$recall.linking)
	closest = min(abs(recall.this-temp$recall.linking))
 	 if (closest<0.02){
    	# identify + filter to row with a recall closest to recall.this
		idx = which.min(abs(recall.this-temp$recall.linking))
    	temp = temp[idx,]
		# concatenate
		if (nrow(df)==0){
			df = temp
		} else {
			df = rbind(df, temp)
		}
	 }
}

# organize for plotting
if (nrow(df)>1) {
	df = dplyr::left_join(df, cp, by="method")
	df$key = paste0(df$pred_name_long, " (", df$Biosample, ")")

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
}

# all pairwise comparisons
if (nrow(df)>1) {
  ## calculate all pair-wise comparison p-values (ref: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1125071/)
  df_p_val = data.frame(t(combn(df$key, 2)))
  colnames(df_p_val) = c("group1", "group2")
  df_p_val$p = 0

  # iterate through rows
  for (i in 1:nrow(df_p_val)) {
    # calculate p-value
    method.1 = dplyr::filter(df, key==df_p_val$group1[i])
    method.2 = dplyr::filter(df, key==df_p_val$group2[i])
    d = log(method.1$enrichment[1]) - log(method.2$enrichment[1]) # difference in log RRs
    SE_d = sqrt(method.1$SE_log_enr[1]**2 + method.2$SE_log_enr[1]**2)
    z = d/SE_d
    p = pnorm(-(abs(z))) * 2 # two-sided p-value
    df_p_val$p[i] = p
  }
  df_p_val$p_adjust = p.adjust(df_p_val$p, method="bonferroni")
  df_p_val$significant = df_p_val$p_adjust < sign_threshold
  print(df_p_val)

  ## plotting
  # format things a little
  df$recall.linking.rounded = round(df$recall.linking, digits=3)
  df$plotting_label = paste0(df$key, " (", df$recall.linking.rounded, ")")
  df = df[order(df$enrichment, decreasing=TRUE),]
  df$plotting_label = factor(df$plotting_label, levels=df$plotting_label, ordered=TRUE)
  df$key = factor(df$key, levels=df$key, ordered=TRUE)

pred_colors = cp$hex
names(pred_colors) = cp$key

  # actually plot
  g = ggplot(data=df, aes(x=key, y=enrichment)) +
    geom_col(aes(fill=key)) +
    geom_errorbar(aes(ymin=CI_enr_low, ymax=CI_enr_high), width=0.2) +
    scale_fill_manual(values=pred_colors, labels=df$plotting_label) +
    labs(fill="Predictor (exact recall)") + ylab(paste0("Enrichment (eQTLs vs. common variants) at recall ", recall.this)) +
    theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), legend.text=element_text(size=7), legend.title=element_text(size=8)) + 
    theme(axis.text.x=element_blank(), axis.title.x=element_blank()) # remove x-axis labels
  


} else {
  df = "Fewer than two predictors achieve this recall."
  df_p_val = "Fewer than two predictors achieve this recall."
  g = ggplot() + theme_void()
}

# save outputs
write.table(df_p_val, out_sign, quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")
write.table(df, out_table, quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")
ggsave(out_plot, g, width=6, height=4)

