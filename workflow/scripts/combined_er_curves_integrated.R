## NEW RECALL
## load libraries
library(dplyr)
library(ggplot2)
library(colorspace)
library(data.table)
library(tidyr)

## INPUTS
methods = snakemake@params$methods %>% strsplit(" ") %>% unlist()
methods_config = fread(snakemake@params$methods_config, sep="\t")
GTExTissue.this = snakemake@wildcards$GTExTissue
cpFilePlotting = snakemake@input$colorPalette
cpList = readRDS(cpFilePlotting)
outDir = snakemake@params$outDir
out_file = snakemake@output$er_combined
outTable_file = snakemake@output$er_combined_table

## gather ER curve tables for methods with that tissue
# if multiple matches, just take first (can adapt to be better later!)
pos_methods = c()
matched_biosamples=c()
print(GTExTissue.this)
for (i in 1:length(methods)){
  methods_temp = dplyr::filter(methods_config, method==methods[i])
  sampleKey_temp = fread(methods_temp$sampleKey[1], sep="\t")
  sampleKey_temp = dplyr::select(sampleKey_temp, biosample, GTExTissue) %>% dplyr::filter(GTExTissue!="")
  if (GTExTissue.this %in% sampleKey_temp$GTExTissue){
    pos_methods = c(pos_methods, methods[i])
    matched_temp = dplyr::filter(sampleKey_temp, GTExTissue==GTExTissue.this)
    matched_biosamples = c(matched_biosamples, matched_temp$biosample[1])
  }
}
print(pos_methods)
print(matched_biosamples)
df_files = data.frame(method=pos_methods, biosample=matched_biosamples)
df_files$file_name = paste0("GTExTissue", GTExTissue.this, ".Biosample", df_files$biosample, ".tsv")
df_files$file_path = file.path(outDir, df_files$method, "enrichmentRecallTables", df_files$file_name)

### FORMAT DATA
## aggregate and process ERCurve tables
# process first table
temp = read.table(file = df_files$file_path[1],header = TRUE, fill = TRUE) 
#temp = drop_na(temp)
method_this = df_files$method[1]
methods_temp = dplyr::filter(methods_config, method==method_this)
temp$method = method_this
temp$pred_name_long = methods_temp$pred_name_long[1]
temp$nPoints = nrow(temp)
ER_curve_all = temp

if (length(pos_methods) > 1) {
  for (i in 2:length(pos_methods)) {
    temp = read.table(file = df_files$file_path[i],header = TRUE, fill = TRUE)
    #temp = drop_na(temp)
    method_this = df_files$method[i]
    methods_temp = dplyr::filter(methods_config, method==method_this)
    temp$method = methods[i]
    temp$pred_name_long = methods_temp$pred_name_long[1]
    temp$nPoints = nrow(temp)
    ER_curve_all = rbind(ER_curve_all, temp)
    }
}

df = ER_curve_all
df[is.na(df)] = 0
df = df[order(df$threshold),]

# data for binary predictors to be plotted as points
df_binary = dplyr::filter(df, nPoints==2, threshold==1)

df = dplyr::filter(df, nPoints!=2)

ylim = 50

g=ggplot(data=df, aes(x=recall.linking, y=enrichment, color=pred_name_long)) +
  geom_line(linewidth=1) +
  geom_point(data=df_binary, aes(x=recall.linking, y=enrichment, color=pred_name_long), size=4) +
  scale_color_manual(values=unlist(cpList)) +
  ylab("Enrichment (eQTLs vs. common variants)") + xlab("Recall (variants overlapping prediction linked to eGene)") +
  labs(col="Predictor") +
  coord_cartesian(ylim=c(0,ylim)) +
  #ylim(c(0, ylim+0.1)) +
  theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8))

ggsave(out_file, g, width=8, height=6)

write.table(df, outTable_file, quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")

