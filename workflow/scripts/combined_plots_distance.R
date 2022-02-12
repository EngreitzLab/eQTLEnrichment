suppressPackageStartupMessages({
  library(plyr)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(egg)
  library(colorspace)
  library(forcats)
})


### INPUT DATA
# methods = c("ABC_20211207_hg38ENCODE_RefSeq_0.01", "ABC_20211207_hg38ENCODE_RefSeq_0.015", "ABC_20211207_hg38ENCODE_RefSeq_0.02",
#             "ChromHMM", "EpiMap",
#             "nearest_gene", "nearest_tss",
#             "dist_to_gene", "dist_to_tss", "within_100kb_of_tss",
#             "reads_by_dist_to_tss", "reads_by_dist_to_tss_norm")
# results_file = "2021_1209_allbiosamples/"
# base = "/oak/stanford/groups/engreitz/Users/sheth/eQTLEnrichment-byBiosample/eQTLEnrichment/workflow/"
# tables = paste0(base, results_file, methods, "/enrichmentTable.tsv")
# varPerTissue = paste0(base, results_file, methods, "/nVariantsPerGTExTissue.tsv")
# predictionMetricsFile_baseline = paste0(base, results_file, "predictionMetrics.tsv") 
# predictionMetricsFile_ABC = paste0(base, results_file, methods[3], "/predictionMetrics.tsv")
# predictionMetricsFile_Other = paste0(base, results_file, "EpiMap", "/predictionMetrics.tsv")

### FOR DISTANCES
method = "ABC_20210712_hg38ENCODE_RefSeq"
distances = c("1000", "10000", "100000", "5000000")
results_file = "2022_0210_test1/"
base = "/oak/stanford/groups/engreitz/Users/sheth/eQTLEnrichment-distance/eQTLEnrichment/workflow/"
tables = paste0(base, results_file, method, "/enrichmentTable.under", distances, "bp.tsv")

### GENERATE PLOTS
## make color palette
cp = data.frame(distance=distances, 
                hex=qualitative_hcl(n=length(distances), palette="Set 2"))
cpList = split(f=cp$distance, x=cp$hex)

## aggregate enrichment tables
enr.all = read.table(tables[1], header = TRUE, fill = TRUE) %>% drop_na()
enr.all$distance = distances[1]

if (length(distances) > 1) {
  for (i in 2:length(distances)) {
    temp = read.table(file = tables[i],header = TRUE, fill = TRUE) %>% drop_na()
    temp$distance = distances[i]
    enr.all = rbind(enr.all, temp)
  }
  print(nrow(enr.all))
}
enr.all = enr.all[order(enr.all$enrichment), ]
max.enr = max(enr.all$enrichment)
enrLabel = 'Enrichment\n(GTEx variants/all common variants)'
enr.all$distance = factor(enr.all$distance, levels=distances, ordered=TRUE)

## enrichment boxplot
boxplot.summary = ggplot(enr.all, aes(x = distance, y = enrichment, fill = distance)) +
  geom_boxplot() +
  theme_minimal() + ylab(enrLabel) + xlab("Variant-gene distance threshold (bp)") +
  scale_fill_manual(values=cpList) +
  ggtitle('Enrichment of eQTLs in predictions') +
  theme(legend.position = 'none') + coord_flip()
  
#pdf(file=paste0(base, results_file, 'boxplot_ABC_bydistance.pdf'), width=7, height=5)
#   print(boxplot.summary)
#dev.off()

## Prediction metrics
pmNames = c("predictionMetrics.under1000bp.tsv", "predictionMetrics.under10000bp.tsv", 
            "predictionMetrics.under100000bp.tsv", "predictionMetrics.under5000000bp.tsv")

predictionMetrics = data.frame('temp','temp','temp', 1, 1) %>% setNames(c("GTExTissue", "metric", "method","value", "distance"))

for (i in 1:length(distances)) {
  df.file = paste0(base, results_file, pmNames[i])
  df = read.table(file=df.file, header=TRUE, sep="\t")
  df$distance = distances[i]
  predictionMetrics = rbind(predictionMetrics, df)
}

predictionMetrics = dplyr::filter(predictionMetrics, GTExTissue!="temp") %>%
  filter(metric %in% c("% variants in any enhancer","Prediction rate given variant in enhancer"))

predictionMetrics$metric = plyr::revalue(predictionMetrics$metric, 
                                    c("% variants in any enhancer"="Overlaps predicted \nenhancer",
                                      "Prediction rate given variant in enhancer"="Overlaps predicted enhancer \nlinked to correct gene"))
small.rates.1 = filter(predictionMetrics, metric=="Overlaps predicted \nenhancer")
small.rates.2 = filter(predictionMetrics, metric=="Overlaps predicted enhancer \nlinked to correct gene")


# ## variant capture and predition rate
# # read in prediction metrics data, set names if necessary
# df1 = read.table(predictionMetricsFile_baseline, header = TRUE, sep = '\t') %>%
#   setNames(c('GTExTissue', 'metric', 'method', 'value'))
# df2 = read.table(predictionMetricsFile_ABC, header = TRUE, sep = '\t') %>%
#   setNames(c('GTExTissue', 'metric', 'method', 'value'))
# df3 = read.table(predictionMetricsFile_Other, header = TRUE, sep = '\t') %>%
#   setNames(c('GTExTissue', 'metric', 'method', 'value'))
# df = rbind(df1, df2, df3)
# df[df == "Cells_EBV-transformed_lymphocytes"] = "LCLs"
# 
# df$method = factor(df$method, levels=methods, ordered=TRUE)
# 
# # subset data from all  metrics
# #small.rates = filter(df, metric %in% c("% variants in any enhancer", "Prediction rate"))
# small.rates = df
# small.rates$metric = plyr::revalue(
#   small.rates$metric,
#   c("% variants in any enhancer" = "Overlaps predicted enhancer",
#     "Prediction rate" = "Overlaps predicted enhancer\nlinked to correct gene",
#     "Prediction rate given variant in enhancer" = "Linked to correct gene\ngiven overlaps predicted enhancer"))
# small.rates.1 = filter(small.rates, metric=="Overlaps predicted enhancer")
# small.rates.2 = filter(small.rates, metric=="Linked to correct gene\ngiven overlaps predicted enhancer")
# 
# small.rates.1$method = fct_rev(small.rates.1$method)
# small.rates.2$method = fct_rev(small.rates.2$method)
# 
# # graph average recall across tissues
# sr.average = ggplot(small.rates, aes(x = metric, y = value)) +
#   geom_bar(stat = "identity", position = "dodge", width = 0.8, aes(fill=method)) +
#   scale_fill_manual(values=cp$hex, name="Method") +
#   theme_minimal() + xlab('') + ylab('Average fraction of fine-mapped, distal noncoding\neQTL variants with PIP >= 0.5') +
#   theme(axis.text.y = element_text(hjust = 1))
# 
# boxplots
sr.overlaps = ggplot(small.rates.1, aes(x = distance, y = value, fill=distance)) +
  geom_boxplot() +
  scale_fill_manual(values=cpList) +
  theme_minimal() + ggtitle('Variants overlapping\npredicted enhancer') +
  ylab('Fraction of fine-mapped, distal noncoding\neQTL variants with PIP >= 0.5') + xlab('') +
  theme(axis.text.y = element_blank(), legend.position='none') + coord_flip()

sr.predicted = ggplot(small.rates.2, aes(x = distance, y = value, fill=distance)) +
  geom_boxplot() +
  scale_fill_manual(values=cpList) +
  theme_minimal() + ggtitle('Variants linked to correct gene\ngiven overlaps predicted enhancer') + xlab('') +
  ylab('Fraction of fine-mapped, distal noncoding\neQTL variants with PIP >= 0.5') +
  theme(axis.text.y = element_blank(), legend.position='none') + coord_flip()
# 
# # horizontal bar plots
# sr.overlaps.bar = ggplot(small.rates.1, aes(x = method, y = value, fill=method)) +
#   geom_bar(stat = "identity", position = "dodge", width = 0.6) +
#   scale_fill_manual(values=cpList) +
#   theme_minimal() + xlab('Overlaps predicted enhancer') + 
#   ylab('Average fraction of fine-mapped, distal\nnoncoding eQTL variants with PIP >= 0.5') +
#   theme(axis.text.y = element_blank(), legend.position='none') + coord_flip()
# 
# sr.predicted.bar = ggplot(small.rates.2, aes(x = method, y = value, fill=method)) +
#   geom_bar(stat = "identity", position = "dodge", width = 0.6) +
#   scale_fill_manual(values=cpList) +
#   theme_minimal() + xlab('Linked to correct gene\ngiven overlaps predicted enhancer') + 
#   ylab('Average fraction of fine-mapped, distal\nnoncoding eQTL variants with PIP >= 0.5') +
#   theme(axis.text.y = element_blank(), legend.position='none') + coord_flip()
# 
# # pdf(file=paste0(base, results_file, 'rates_overlap.pdf'), width=7, height=5)
# #   print(sr.overlaps)
# # dev.off()
# # 
# # pdf(file=paste0(base, results_file, 'rates_predicted.pdf'), width=7, height=5)
# # print(sr.predicted)
# # dev.off()
# 
# # combine and save
# final = ggarrange(boxplot.summary, sr.average, nrow=1, ncol=2)
# # pdf(file=paste0(base, results_file, 'summary_figure_combined_v2.pdf'), width=12, height=5)
# #   print(final)
# # dev.off()
# 
# final.v2 =  ggarrange(boxplot.v2, sr.overlaps, sr.predicted, nrow=1, ncol=3)
#  pdf(file=paste0(base, results_file, 'combined_figure_horizontal.v2.pdf'), width=12, height=5)
#     print(final.v2)
#  dev.off()
# 
# final.v3 =  ggarrange(boxplot.v2, sr.overlaps.bar, sr.predicted.bar, nrow=1, ncol=3)
# #pdf(file=paste0(base, results_file, 'combined_figure_horizontal_bars.pdf'), width=14, height=5)
# #    print(final.v3)
# #dev.off()

test.final = ggarrange(boxplot.summary, sr.overlaps, sr.predicted, nrow=1, ncol=3)
 pdf(file=paste0(base, results_file, 'combined_figure_by_distance.pdf'), width=12, height=5)
    print(test.final)
 dev.off()