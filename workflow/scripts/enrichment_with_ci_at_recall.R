## NEW RECALL
## load libraries
library(dplyr)
library(ggplot2)
library(colorspace)
library(data.table)
library(tidyr)
library(ggpubr)

## INPUTS
methods = snakemake@params$methods %>% strsplit(" ") %>% unlist()
methods_config = fread(snakemake@params$methods_config, sep="\t")
GTExTissue.this = snakemake@wildcards$GTExTissue
recall.this = snakemake@wildcards$recall %>% as.numeric()
cpFilePlotting = snakemake@input$colorPalette
cpList = readRDS(cpFilePlotting)
outDir = snakemake@params$outDir
out_file = snakemake@output$enr_at_recall
outTable_file = snakemake@output$enr_at_recall_table
outSign_file = snakemake@output$sign_table

sign_threshold = 0.05

## gather ER curve tables for methods with this GTEx tissue
# if multiple biosample matches, just take first (can adapt to be better later!)
pos_methods = c()
matched_biosamples=c()
for (i in 1:length(methods)){
  methods_temp = dplyr::filter(methods_config, method==methods[i])
  sampleKey_temp = fread(methods_temp$sampleKey[1], sep="\t")
  sampleKey_temp = dplyr::select(sampleKey_temp, biosample, GTExTissue) %>% dplyr::filter(GTExTissue!="")
  if (GTExTissue.this %in% sampleKey_temp$GTExTissue) {
    pos_methods = c(pos_methods, methods[i])
    matched_temp = dplyr::filter(sampleKey_temp, GTExTissue==GTExTissue.this)
    matched_biosamples = c(matched_biosamples, matched_temp$biosample[1])
  }
}

df_files = data.frame(method=pos_methods, biosample=matched_biosamples)
df_files$file_name = paste0("GTExTissue", GTExTissue.this, ".Biosample", df_files$biosample, ".tsv")
df_files$file_path = file.path(outDir, df_files$method, "enrichmentRecallTables", df_files$file_name)

## read in ER tables and pull enrichment, recall, and CI bounds for threshold yielding recall closest to recall.this
for (i in 1:length(pos_methods)) {
  temp = read.table(file = df_files$file_path[i], header=TRUE, fill=TRUE)
  temp = drop_na(temp)
  method_this = df_files$method[i] # could also do pos_methods[i]?
  methods_temp = dplyr::filter(methods_config, method==method_this)
  temp$method = method_this
  temp$pred_name_long = methods_temp$pred_name_long[1]
  
  # does the range of recalls overlap with recall.this? only continue if so (this means also filtering out binary predictors)
  n_thresholds = nrow(temp)
  recall_low = min(temp$recall.linking)
  recall_high = max(temp$recall.linking)
  if (n_thresholds>2 & (recall.this>recall_low & recall.this<recall_high)){
    # identify + filter to row with a recall closest to recall.this
    recall_closest = temp$recall.linking[which.min(abs(recall.this-temp$recall.linking))]
    temp = dplyr::filter(temp, recall.linking==recall_closest)

    # concatenate
    if (exists(quote(df_full))) {
      df_full = rbind(df_full, temp)
    } else {
      df_full = temp
    }
  }
}

if (exists(quote(df_full))) {
  ## calculate all pair-wise comparison p-values (ref: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1125071/)
  df_p_val = data.frame(t(combn(df_full$pred_name_long, 2)))
  colnames(df_p_val) = c("group1", "group2")
  df_p_val$p = 0
  df_p_val$y.position = 0

  # iterate through rows
  for (i in 1:nrow(df_p_val)) {
    # calculate p-value
    method.1 = dplyr::filter(df_full, pred_name_long==df_p_val$group1[i])
    method.2 = dplyr::filter(df_full, pred_name_long==df_p_val$group2[i])
    d = log(method.1$enrichment[1]) - log(method.2$enrichment[1]) # difference in log RRs
    SE_d = sqrt(method.1$SE_log_enr[1]**2 + method.2$SE_log_enr[1]**2)
    z = d/SE_d
    p = pnorm(-(abs(z))) * 2 # two-sided p-value
    df_p_val$p[i] = p
  }
  print(df_p_val)

  # filter to sign threshold for plotting
  df_p_val_sig = dplyr::filter(df_p_val, p < sign_threshold)
  print(df_p_val_sig)

  # calculate y-positions of brackets
  if (nrow(df_p_val_sig)> 0) {
    for (i in 1:nrow(df_p_val_sig)) {
      method.1 = dplyr::filter(df_full, pred_name_long==df_p_val_sig$group1[i])
      method.2 = dplyr::filter(df_full, pred_name_long==df_p_val_sig$group2[i])
      y = round(max(method.1$CI_enr_high[1], method.2$CI_enr_high[1]) + 3, digits=0)
      print(y)
      while (y %in% round(df_p_val_sig$y.position, digits=0)) {
        y = y + 1
      }
      df_p_val_sig$y.position[i] = y
    }
  }

  print(df_p_val_sig)

  ## plotting
  # format things a little
  df_full$recall.linking.rounded = round(df_full$recall.linking, digits=3)
  df_full$plotting_label = paste0(df_full$pred_name_long, " (", df_full$recall.linking.rounded, ")")
  df_full = df_full[order(df_full$enrichment, decreasing=TRUE),]
  df_full$plotting_label = factor(df_full$plotting_label, levels=df_full$plotting_label, ordered=TRUE)
  df_full$pred_name_long = factor(df_full$pred_name_long, levels=df_full$pred_name_long, ordered=TRUE)

  # actually plot
  g = ggplot(data=df_full, aes(x=pred_name_long, y=enrichment)) +
    geom_col(aes(fill=pred_name_long)) +
    geom_errorbar(aes(ymin=CI_enr_low, ymax=CI_enr_high), width=0.2) +
    scale_fill_manual(values=unlist(cpList), labels=df_full$plotting_label) +
    labs(fill="Predictor (exact recall)") + ylab(paste0("Enrichment (eQTLs vs. common variants) at recall of ", recall.this)) +
    theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8)) + 
    theme(axis.text.x=element_blank(), axis.title.x=element_blank()) # remove x-axis labels
  


} else {
  df_full = " "
  df_p_val = " "
  g = ggplot() + theme_void()
}

# save outputs
write.table(df_p_val, outSign_file, quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")
write.table(df_full, outTable_file, quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")
ggsave(out_file, g, width=6, height=4)

