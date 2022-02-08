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
methods = c("ABC_20211207_hg38ENCODE_RefSeq_0.01", "ABC_20211207_hg38ENCODE_RefSeq_0.015", "ABC_20211207_hg38ENCODE_RefSeq_0.02",
            "ChromHMM", "EpiMap",
            "nearest_gene", "nearest_tss",
            "dist_to_gene", "dist_to_tss", "within_100kb_of_tss",
            "reads_by_dist_to_tss", "reads_by_dist_to_tss_norm")
results_file = "2021_1209_allbiosamples/"
base = "/oak/stanford/groups/engreitz/Users/sheth/eQTLEnrichment-byBiosample/eQTLEnrichment/workflow/"
tables = paste0(base, results_file, methods, "/enrichmentTable.tsv")
varPerTissue = paste0(base, results_file, methods, "/nVariantsPerGTExTissue.tsv")
predictionMetricsFile_baseline = paste0(base, results_file, "predictionMetrics.tsv") 
predictionMetricsFile_ABC = paste0(base, results_file, methods[3], "/predictionMetrics.tsv")
predictionMetricsFile_Other = paste0(base, results_file, "EpiMap", "/predictionMetrics.tsv")

### Tissue/biosample keys
ABC.key = data.frame(GTExTissue = c("Brain_Cortex", "Prostate", "Muscle_Skeletal", "Artery_Tibial", "Skin_Not_Sun_Exposed_Suprapubic", "Brain_Cerebellum", "Pancreas", "Spleen", "Adrenal_Gland", "Testis", "Lung", "Skin_Sun_Exposed_Lower_leg", "Heart_Atrial_Appendage", "Thyroid", "Colon_Sigmoid", "Breast_Mammary_Tissue", "Artery_Aorta", "Cells_EBV-transformed_lymphocytes", "Heart_Left_Ventricle", "Artery_Coronary", "Nerve_Tibial", "Liver", "Whole_Blood", "Colon_Transverse", "Stomach"),
                     biosample = c("middle_frontal_area_46_ID_414", "epithelial_cell_of_prostate_ID_2659", "muscle_of_trunk_ID_2601", "tibial_artery_ID_2651", "suprapubic_skin_ID_2750", "neuronal_stem_cell_ID_2716", "pancreas_ID_274", "spleen_ID_0", "adrenal_gland_ID_2188", "testis_ID_2718", "lung_ID_2152", "keratinocyte_ID_2339", "right_atrium_auricular_region_ID_2728", "thyroid_gland_ID_2682", "sigmoid_colon_ID_2752", "mammary_epithelial_cell_ID_2625", "aorta_ID_2340", "GM12878_ID_2724", "heart_left_ventricle_ID_169", "coronary_artery_ID_2664", "tibial_nerve_ID_2180", "liver_ID_411", "CD4-positive_alpha-beta_T_cell_ID_119", "transverse_colon_ID_2733", "stomach_ID_2411"))
ABC.key$combined = paste0(ABC.key$GTExTissue, ".", ABC.key$biosample)
ChromHMM.key = data.frame(GTExTissue = c("Brain_Cortex", "Muscle_Skeletal", "Skin_Not_Sun_Exposed_Suprapubic", "Brain_Cerebellum", "Pancreas", "Spleen", "Adrenal_Gland", "Lung", "Skin_Sun_Exposed_Lower_leg", "Heart_Atrial_Appendage", "Colon_Sigmoid", "Breast_Mammary_Tissue", "Artery_Aorta", "Cells_EBV-transformed_lymphocytes", "Heart_Left_Ventricle", "Liver", "Whole_Blood", "Stomach", "Esophagus_Muscularis", "Adipose_Subcutaneous", "Brain_Nucleus_accumbens_basal_ganglia"), 
                          biosample=c("E073", "E108", "E056", "E007", "E098", "E113", "E080", "E096", "E126", "E104", "E106", "E119", "E065", "E116", "E095", "E066", "E034", "E092", "E079", "E023", "E074"))
ChromHMM.key$combined = paste0(ChromHMM.key$GTExTissue, ".", ChromHMM.key$biosample)
EpiMap.key = data.frame(GTExTissue = c("Brain_Cortex", "Prostate", "Muscle_Skeletal", "Artery_Tibial", "Skin_Not_Sun_Exposed_Suprapubic", "Brain_Cerebellum", "Pancreas", "Spleen", "Adrenal_Gland", "Testis", "Lung", "Skin_Sun_Exposed_Lower_leg", "Heart_Atrial_Appendage", "Thyroid", "Colon_Sigmoid", "Breast_Mammary_Tissue", "Artery_Aorta", "Cells_EBV-transformed_lymphocytes", "Heart_Left_Ventricle", "Artery_Coronary", "Nerve_Tibial", "Liver", "Whole_Blood", "Colon_Transverse", "Stomach", "Esophagus_Muscularis", "Adipose_Subcutaneous", "Esophagus_Mucosa", "Esophagus_Gastroesophageal_Junction"),	
                        biosample=c("BSS01271", "BSS01459", "BSS01572", "BSS01839", "BSS01678", "BSS01372", "BSS01408", "BSS01632", "BSS00054", "BSS01715", "BSS01190", "BSS01180", "BSS01506", "BSS01831", "BSS01542", "BSS01209", "BSS00079", "BSS00439", "BSS00512", "BSS00242", "BSS01842", "BSS00511", "BSS00190", "BSS01848", "BSS01639", "BSS00316", "BSS01671", "BSS00323", "BSS00382"))
EpiMap.key$combined = paste0(EpiMap.key$GTExTissue, ".", EpiMap.key$biosample)

### GENERATE PLOTS
## make color palette
cp = data.frame(method=methods, 
                hex=qualitative_hcl(n=length(methods), palette="Set 2"))
cpList = split(f=cp$method, x=cp$hex)

## aggregate enrichment tables
enr.all = read.table(tables[1], header = TRUE, fill = TRUE) %>% drop_na()
enr.all$method = methods[1]
enr.all$tissue.biosample = paste0(enr.all$GTExTissue, ".", enr.all$Biosample)
enr.all = filter(enr.all, tissue.biosample %in% ABC.key$combined)
print(nrow(enr.all))

if (length(methods) > 1) {
  for (i in 2:length(methods)) {
    temp = read.table(file = tables[i],header = TRUE, fill = TRUE) %>% drop_na()
    temp$method = methods[i]
    temp$tissue.biosample = paste0(temp$GTExTissue, ".", temp$Biosample)
    if (methods[i] == "ChromHMM"){
      key = ChromHMM.key
    } else if (methods[i] == "EpiMap"){
      key = EpiMap.key
    } else {
      key = ABC.key
    }
    temp = filter(temp, tissue.biosample %in% key$combined)
    enr.all = rbind(enr.all, temp)
  }
  print(nrow(enr.all))
}
enr.all = enr.all[order(enr.all$enrichment), ]
max.enr = max(enr.all$enrichment)
enrLabel = 'Enrichment\n(GTEx variants/all common variants)'
enr.all$method = factor(enr.all$method, levels=methods, ordered=TRUE)

## enrichment boxplot
boxplot.summary = ggplot(enr.all, aes(x = method, y = enrichment, fill = method)) +
  geom_boxplot() +
  theme_minimal() + ylab(enrLabel) + xlab('') +
  scale_fill_manual(values=cpList) +
  theme(legend.position = 'none', axis.text.x=element_blank())

enr.all$method = fct_rev(enr.all$method)
boxplot.v2 = ggplot(enr.all, aes(x = method, y = enrichment, fill = method)) +
  geom_boxplot() +
  coord_flip() +
  theme_minimal() + ylab(enrLabel) + xlab('') +
  scale_fill_manual(values=cpList) +
  theme(legend.position = 'none') +
  ggtitle("Enrichment of eQTLs in predictions\n for matching tissues & biosamples")
  
pdf(file=paste0(base, results_file, 'boxplot_horizontal.pdf'), width=7, height=5)
   print(boxplot.v2)
dev.off()

## variant capture and predition rate
# read in prediction metrics data, set names if necessary
df1 = read.table(predictionMetricsFile_baseline, header = TRUE, sep = '\t') %>%
  setNames(c('GTExTissue', 'metric', 'method', 'value'))
df2 = read.table(predictionMetricsFile_ABC, header = TRUE, sep = '\t') %>%
  setNames(c('GTExTissue', 'metric', 'method', 'value'))
df3 = read.table(predictionMetricsFile_Other, header = TRUE, sep = '\t') %>%
  setNames(c('GTExTissue', 'metric', 'method', 'value'))
df = rbind(df1, df2, df3)
df[df == "Cells_EBV-transformed_lymphocytes"] = "LCLs"

df$method = factor(df$method, levels=methods, ordered=TRUE)

# subset data from all  metrics
#small.rates = filter(df, metric %in% c("% variants in any enhancer", "Prediction rate"))
small.rates = df
small.rates$metric = plyr::revalue(
  small.rates$metric,
  c("% variants in any enhancer" = "Overlaps predicted enhancer",
    "Prediction rate" = "Overlaps predicted enhancer\nlinked to correct gene",
    "Prediction rate given variant in enhancer" = "Linked to correct gene\ngiven overlaps predicted enhancer"))
small.rates.1 = filter(small.rates, metric=="Overlaps predicted enhancer")
small.rates.2 = filter(small.rates, metric=="Linked to correct gene\ngiven overlaps predicted enhancer")

small.rates.1$method = fct_rev(small.rates.1$method)
small.rates.2$method = fct_rev(small.rates.2$method)

# graph average recall across tissues
sr.average = ggplot(small.rates, aes(x = metric, y = value)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.8, aes(fill=method)) +
  scale_fill_manual(values=cp$hex, name="Method") +
  theme_minimal() + xlab('') + ylab('Average fraction of fine-mapped, distal noncoding\neQTL variants with PIP >= 0.5') +
  theme(axis.text.y = element_text(hjust = 1))

# boxplots
sr.overlaps = ggplot(small.rates.1, aes(x = method, y = value, fill=method)) +
  geom_boxplot() +
  scale_fill_manual(values=cpList) +
  theme_minimal() + ggtitle('Variants overlapping\npredicted enhancer') + 
  ylab('Fraction of fine-mapped, distal noncoding\neQTL variants with PIP >= 0.5') + xlab('') +
  theme(axis.text.y = element_blank(), legend.position='none') + coord_flip()

sr.predicted = ggplot(small.rates.2, aes(x = method, y = value, fill=method)) +
  geom_boxplot() +
  scale_fill_manual(values=cpList) +
  theme_minimal() + ggtitle('Variants linked to correct gene\ngiven overlaps predicted enhancer') + xlab('') +
  ylab('Fraction of fine-mapped, distal noncoding\neQTL variants with PIP >= 0.5') +
  theme(axis.text.y = element_blank(), legend.position='none') + coord_flip()

# horizontal bar plots
sr.overlaps.bar = ggplot(small.rates.1, aes(x = method, y = value, fill=method)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  scale_fill_manual(values=cpList) +
  theme_minimal() + xlab('Overlaps predicted enhancer') + 
  ylab('Average fraction of fine-mapped, distal\nnoncoding eQTL variants with PIP >= 0.5') +
  theme(axis.text.y = element_blank(), legend.position='none') + coord_flip()

sr.predicted.bar = ggplot(small.rates.2, aes(x = method, y = value, fill=method)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  scale_fill_manual(values=cpList) +
  theme_minimal() + xlab('Linked to correct gene\ngiven overlaps predicted enhancer') + 
  ylab('Average fraction of fine-mapped, distal\nnoncoding eQTL variants with PIP >= 0.5') +
  theme(axis.text.y = element_blank(), legend.position='none') + coord_flip()

# pdf(file=paste0(base, results_file, 'rates_overlap.pdf'), width=7, height=5)
#   print(sr.overlaps)
# dev.off()
# 
# pdf(file=paste0(base, results_file, 'rates_predicted.pdf'), width=7, height=5)
# print(sr.predicted)
# dev.off()

# combine and save
final = ggarrange(boxplot.summary, sr.average, nrow=1, ncol=2)
# pdf(file=paste0(base, results_file, 'summary_figure_combined_v2.pdf'), width=12, height=5)
#   print(final)
# dev.off()

final.v2 =  ggarrange(boxplot.v2, sr.overlaps, sr.predicted, nrow=1, ncol=3)
 pdf(file=paste0(base, results_file, 'combined_figure_horizontal.v2.pdf'), width=12, height=5)
    print(final.v2)
 dev.off()

final.v3 =  ggarrange(boxplot.v2, sr.overlaps.bar, sr.predicted.bar, nrow=1, ncol=3)
#pdf(file=paste0(base, results_file, 'combined_figure_horizontal_bars.pdf'), width=14, height=5)
#    print(final.v3)
#dev.off()

