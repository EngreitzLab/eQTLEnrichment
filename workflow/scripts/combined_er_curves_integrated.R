## NEW RECALL
## load libraries
library(dplyr)
library(ggplot2)
library(colorspace)
library(data.table)
library(tidyr)

## INPUTS
methods = snakemake@params$methods %>% strsplit(" ") %>% unlist()
GTExTissue.this = snakemake@wildcards$GTExTissue
cpFilePlotting = snakemake@input$colorPalette
cpList = readRDS(cpFilePlotting)
methods_config = read.table(file=snakemake@params$methods_config, header=TRUE, sep="\t")
outDir = snakemake@params$outDir
out_file = snakemake@output$er_combined
outTable_file = snakemake@output$er_combined_table

## gather ER curve tables for methods with that tissue
# format methods_config columns
methods_config$GTExTissue_map = substr(methods_config$GTExTissue_map, 3, (nchar(methods_config$GTExTissue_map)-2))
methods_config$GTExTissue_map = gsub('"', "", methods_config$GTExTissue_map)
methods_config$GTExTissue_map = gsub(',', "", methods_config$GTExTissue_map)
methods_config$GTExTissue_map = strsplit(methods_config$GTExTissue_map, " ")

methods_config$biosample_map = substr(methods_config$biosample_map, 3, (nchar(methods_config$biosample_map)-2))
methods_config$biosample_map = gsub('"', "", methods_config$biosample_map)
methods_config$biosample_map = gsub(',', "", methods_config$biosample_map)
methods_config$biosample_map = strsplit(methods_config$biosample_map, " ")

# if multiple matches, just take first (can adapt to be better later!)
pos_methods = c()
matched_biosamples=c()
for (i in length(methods)){
  methods_temp = dplyr::filter(methods_config, method==methods[i])
  
  GTExTissues_temp = methods_temp$GTExTissue_map[i]
  if (GTExTissue.this %in% GTExTissues_temp){
    pos_methods = c(pos_methods, methods[i])
    
    matches_temp = data.frame(GTExTissue=methods_temp$GTExTissue_map[i], Biosample=methods_temp$biosample_map[i])
    matches_temp = dplyr::filter(matches_temp, GTExTissue==GTExTissue.this)
    matched_biosamples = c(matched_biosamples, matched_temp$Biosample[1])
  }
}

df_files = data.frame(method=post_methods, biosample=matched_biosamples)
df_files$file_name = paste0("GTExTissue", GTExTissue.this, ".Biosample", df_files$biosample, ".tsv")
df_files$file_path = file.path(outDir, "EnrichmentRecallTables", df_files$file_name)

### FORMAT DATA
## aggregate and process ERCurve tables
# process first table
temp = read.table(file = df_files$file_path[1],header = TRUE, fill = TRUE) %>% drop_na()
method_this = df_files$method[1]
methods_temp = dplyr::filter(methods_config, method==method_this)
temp$method = method_this
temp$pred_name_long = methods_temp$pred_name_long[1]
temp$nPoints = nrow(temp)
ER_curve_all = temp

if (length(methods) > 1) {
  for (i in 2:length(methods)) {
    temp = read.table(file = df_files$file_path[i],header = TRUE, fill = TRUE) %>% drop_na()
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

df_binary = dplyr::filter(df, nPoints==2, threshold==1)

df = dplyr::filter(df, nPoints!=2)

# data for binary predictors to be plotted as points

g=ggplot(data=df, aes(x=recall.LCL, y=enrichment, color=pred_name_long)) +
  geom_line(linewidth=1) +
  geom_point(data=df_binary, aes(x=recall.linking, y=enrichment, color=pred_name_long), size=4) +
  scale_color_manual(values=unlist(cpList)) +
  ylab("Enrichment (eQTLs vs. common variants)") + xlab("Recall (variants overlapping prediction linked to eGene") +
  labs(col="Predictor") +
  ylim(c(0, 50)) +
  theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8))

ggsave(out_file, g, width=8, height=6)

write.table(df, outTable_file, quote=FALSE, col.names=TRUE, sep="\t")

