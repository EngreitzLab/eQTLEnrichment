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
distances = snakemake@params$distances %>% as.character %>% strsplit(" ") %>% unlist()
enrTableFiles = snakemake@input$enrichmentTables %>% strsplit(" ") %>% unlist()
methods_config = fread(file=snakemake@params$methods_config, sep="\t")

cpFilePlotting = snakemake@input$colorPalette
cpList = readRDS(cpFilePlotting)

outFile = snakemake@output$outFile
outEnrAll = snakemake@output$enrAllTable
outPredMetrics = snakemake@output$predictionMetrics
outDir = snakemake@params$outDir

### FORMAT DATA
## aggregate and process enrichment tables
# first one
temp = read.table(file = enrTableFiles[1],header = TRUE, fill = TRUE)
temp$method = methods[1]
temp$distance = distances[1]
temp$tissue.biosample = paste0(temp$GTExTissue, ".", temp$Biosample)
methods_temp = dplyr::filter(methods_config, method==methods[1])
sampleKey_temp = fread(methods_temp$sampleKey[1], sep="\t")
sampleKey_temp =  dplyr::select(sampleKey_temp, biosample, GTExTissue) %>% dplyr::filter(GTExTissue!="")
sampleKey_temp$key = paste0(sampleKey_temp$GTExTissue, ".", sampleKey_temp$biosample)

temp$pred_name_long = methods_temp$pred_name_long[1]
temp = filter(temp, tissue.biosample %in% sampleKey_temp$key)
enr.all = temp

# loop through rest, filtering to tissue/biosample matches
  for (i in 1:length(methods)) {
    for (j in 1:length(distances)){
    temp = read.table(file = enrTableFiles[1+(i-1)*length(distances)+j-1],header = TRUE, fill = TRUE)
    temp$method = methods[i]
    temp$distance = distances[j]
    temp$tissue.biosample = paste0(temp$GTExTissue, ".", temp$Biosample)
    methods_temp = dplyr::filter(methods_config, method==methods[i])
    sampleKey_temp = fread(methods_temp$sampleKey[1], sep="\t")
    sampleKey_temp =  dplyr::select(sampleKey_temp, biosample, GTExTissue) %>% dplyr::filter(GTExTissue!="")
    sampleKey_temp$key = paste0(sampleKey_temp$GTExTissue, ".", sampleKey_temp$biosample)
    temp$pred_name_long = methods_temp$pred_name_long[1]
    temp = filter(temp, tissue.biosample %in% sampleKey_temp$key)
    enr.all = rbind(enr.all, temp)
    }
  }
enr.all = distinct(enr.all)

# get stats, edit names, define labels 
enr.all = enr.all[order(enr.all$enrichment), ]
max.enr = max(enr.all$enrichment[is.finite(enr.all$enrichment)])
enrLabel = 'Enrichment\n(GTEx variants/all common variants)'

enr.all$method = factor(enr.all$method, levels=methods, ordered=TRUE)

## aggregate and process prediction metrics
# predTable = os.path.join(config["outDir"], "{method}", "predictionTables", "GTExTissue{GTExTissue}.Biosample{Biosample}.byDistance.tsv")
# first one
method.this = methods[1]
methods_temp = dplyr::filter(methods_config, method==method.this)

sampleKey_temp = fread(methods_temp$sampleKey[1], sep="\t")
sampleKey_temp =  dplyr::select(sampleKey_temp, biosample, GTExTissue) %>% dplyr::filter(GTExTissue!="")
sampleKey_temp$key = paste0(sampleKey_temp$GTExTissue, ".", sampleKey_temp$biosample)
key_file = paste0("GTExTissue", sampleKey_temp$GTExTissue[1], ".", "Biosample", sampleKey_temp$biosample[1], ".byDistance.tsv")
key_path = file.path(outDir, method.this, "predictionTables", key_file)
temp = read.table(file = key_path, header = TRUE, fill = TRUE) %>% drop_na()
temp$GTExTissue = sampleKey_temp$GTExTissue[1]
temp$Biosample = sampleKey_temp$biosample[1]
temp$method = method.this
pred.all = temp

for (i in 1:length(methods)){
  method.this = methods[i]
  methods_temp = dplyr::filter(methods_config, method==method.this)
  sampleKey_temp = fread(methods_temp$sampleKey[1], sep="\t")
  sampleKey_temp =  dplyr::select(sampleKey_temp, biosample, GTExTissue) %>% dplyr::filter(GTExTissue!="")
  sampleKey_temp$key = paste0(sampleKey_temp$GTExTissue, ".", sampleKey_temp$biosample)
  print(sampleKey_temp)
    for (j in 1:nrow(sampleKey_temp)){
      print(sampleKey_temp)
      print(j)
      print(sampleKey_temp$GTExTissue[j])
      print(sampleKey_temp$biosample[j])
      key_file = paste0("GTExTissue", sampleKey_temp$GTExTissue[j], ".", "Biosample", sampleKey_temp$biosample[j], ".byDistance.tsv")
      key_path = file.path(outDir, method.this, "predictionTables", key_file)
      temp = read.table(file = key_path, header = TRUE, fill = TRUE) %>% drop_na()
      temp$GTExTissue = sampleKey_temp$GTExTissue[j]
      temp$Biosample = sampleKey_temp$biosample[j]
      temp$method = method.this
      pred.all = rbind(pred.all, temp)
    }
}
pred.all = distinct(pred.all)

# add column with plotting name by left_joining with methods_config selecting just the method and plotting name columns
methods_config.select = dplyr::select(methods_config, method, pred_name_long)
pred.all = left_join(pred.all, methods_config.select, by=c("method"))
#enr.all = left_join(enr.all, methods_config.select, by=c("method"))

# make y-axis labels (distance range, min-max variants)
pred.all$distance.label = 0
enr.all$distance.label = 0

distances_plus0 = c(0, distances)
for (i in 2:length(distances_plus0)){
  distance.this = distances_plus0[i]
  if (distance.this == 30000000){
    label = "All variants\n"
  } else {
    label = paste0(distances_plus0[i-1], "-", distance.this, "bp\n")
  }
  pred.filt = dplyr::filter(pred.all, distance==distance.this)
  label.full = paste0(label, "N = ", min(pred.filt$total.variants), "-", max(pred.filt$total.variants), " variants")
  pred.all$distance.label[pred.all$distance==distance.this] = label.full
  enr.all$distance.label[enr.all$distance==distance.this] = label.full
}

# for biosamples_predictors plot 
# ordered.methods = c('In element (DHS) & closest gene', 'ABC_A=DNase, C=Avg. Intact Hi-C', 'ENCODE-E2G', 'EpiMap', 'EPIraction', 'ABC_A=DNase x H3K27ac, C=Avg. Intact Hi-C')
# enr.all$pred_name_long = factor(enr.all$pred_name_long, levels=ordered.methods)
# pred.all$pred_name_long = factor(pred.all$pred_name_long, levels=ordered.methods)

### GENERATE PLOTS
## enrichment
enr.boxplot = ggplot(enr.all, aes(x = distance.label, y = enrichment, fill = pred_name_long)) +
  geom_boxplot() +
  coord_flip() +
  theme_minimal() + ylab("Enrichment\n(GTEx variants/all common variants)") + xlab('') +
  scale_fill_manual(values=cpList) +
  theme(legend.position = 'none') +
  ggtitle("Enrichment of variants\nin predicted enhancers")

## overlaps predicted enhancer : recall.total
sr.overlaps = ggplot(pred.all, aes(x = distance.label, y = recall.total, fill=pred_name_long)) +
  geom_boxplot() +
  scale_fill_manual(values=cpList) +
  theme_minimal() + ggtitle('Variants overlapping\npredicted enhancers') + 
  ylab('Fraction of GTEx variants') + xlab('') +
  theme(axis.text.y = element_blank(), legend.position='none') + coord_flip()

## linked to correct eGene
sr.predicted = ggplot(pred.all, aes(x = distance.label, y = correctGene.ifOverlap, fill=pred_name_long)) +
  geom_boxplot() +
  scale_fill_manual(values=cpList) +
  theme_minimal() + ggtitle('Variants linked to correct gene,\ngiven overlapping predicted enhancer') + xlab('') +
  ylab('Fraction of GTEx variants\noverlapping predicted enhancers') +
  theme(axis.text.y = element_blank()) + coord_flip()

## save final plots
pdf(file=outFile, width=12, height=5)
  all.tissues =  ggarrange(enr.boxplot, sr.overlaps, sr.predicted, nrow=1, ncol=3)
dev.off()

pred.all = dplyr::select(pred.all, -distance.label)
enr.all = dplyr::select(enr.all, -distance.label)

write.table(pred.all, outPredMetrics, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
write.table(enr.all, outEnrAll, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
