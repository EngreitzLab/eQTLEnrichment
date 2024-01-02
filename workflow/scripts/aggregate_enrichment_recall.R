suppressPackageStartupMessages({
  library(plyr)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(egg)
  library(colorspace)
  library(forcats)
  library(data.table)
})


### INPUTS
methods = snakemake@params$methods %>% strsplit(" ") %>% unlist()
methods_config = fread(file=snakemake@params$methods_config, sep="\t")
GTExTissues = snakemake@params$matchedGTExTissues %>% strsplit(" ") %>% unlist()

cpFilePlotting = snakemake@input$colorPalette
cpList = readRDS(cpFilePlotting)
outDir = snakemake@params$outDir

outFilePlot = snakemake@output$er_aggregate_plot
outFileTable = snakemake@output$er_aggregate_table

### FORMAT DATA-- goal: table with (for each method): threshold, aggregate recall, aggregate enrichment
# aggregate and process enrichment-recall tables (per tissue)
for (i in 1:length(GTExTissues)) {
  GTExTissue.this = GTExTissues[i]
  fname = paste0("enrichmentRecall.GTExTissue", GTExTissue.this, ".tsv")
  fpath = file.path(outDir, "plots", fname)
  print(fpath)
  temp = read.table(fpath, header=TRUE, sep="\t")
  temp$GTExTissue = GTExTissue.this
  temp$variants.linking = temp$total.variants * temp$recall.linking

  if (i==1){
    df = temp
  } else {
    df = rbind(df, temp)
  }
}

# group by method(pred_name_long)/threshold: enrichment=mean; correct.var=sum, total.variants=sum
df_agg = df %>% group_by(pred_name_long, threshold) %>%
  summarize(enrichment = mean(enrichment), variants.linking = sum(variants.linking), total.variants = sum(total.variants), nPoints = mean(nPoints))

df_agg$recall.linking = df_agg$variants.linking/df_agg$total.variants

### PLOTTING
# data for binary predictors to be plotted as points
df_binary = dplyr::filter(df_agg, nPoints==2, threshold==1)

df_cont = dplyr::filter(df_agg, nPoints!=2)

ylim = 50

g=ggplot(data=df_cont, aes(x=recall.linking, y=enrichment, color=pred_name_long)) +
  geom_line(linewidth=1) +
  geom_point(data=df_binary, aes(x=recall.linking, y=enrichment, color=pred_name_long), size=4) +
  scale_color_manual(values=unlist(cpList)) +
  ylab("Enrichment (eQTLs vs. common variants)") + xlab("Recall (variants overlapping prediction linked to eGene)") +
  labs(col="Predictor") +
  coord_cartesian(ylim=c(0,ylim)) +
  #ylim(c(0, ylim+0.1)) +
  theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8))

ggsave(outFilePlot, g, width=8, height=6)

write.table(df_agg, outFileTable, quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")
