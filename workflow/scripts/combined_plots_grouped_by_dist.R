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

# access via config directly? add methods_config to input
methods_config = fread(file=snakemake@input$methods_config, header=TRUE, sep="\t")
methods_config$GTExTissue_map = substr(methods_config$GTExTissue_map, 3, (nchar(methods_config$GTExTissue_map)-2))
methods_config$GTExTissue_map = gsub('"', "", methods_config$GTExTissue_map)
methods_config$GTExTissue_map = gsub(',', "", methods_config$GTExTissue_map)
methods_config$GTExTissue_map = strsplit(methods_config$GTExTissue_map, " ")

methods_config$biosample_map = substr(methods_config$biosample_map, 3, (nchar(methods_config$biosample_map)-2))
methods_config$biosample_map = gsub('"', "", methods_config$biosample_map)
methods_config$biosample_map = gsub(',', "", methods_config$biosample_map)
methods_config$biosample_map = strsplit(methods_config$biosample_map, " ")

cpFilePlotting = snakemake@input$colorPalette
cpList = readRDS(cpFilePlotting)

outFile = snakemake@output$er_combined
outEnrAll = snakemake@output$enrAllTable
outPredMetrics = snakemake@output$predictionMetrics
outDir = snakemake@params$outDir

### FORMAT DATA
## aggregate and process enrichment tables
# first one
temp = read.table(file = enrTableFiles[1],header = TRUE, fill = TRUE) %>% drop_na()
temp$method = methods[1]
temp$distance = distances[1]
temp$tissue.biosample = paste0(temp$GTExTissue, ".", temp$Biosample)
methods_temp = dplyr::filter(methods_config, method==methods[1])
GTEx_key = methods_temp$GTExTissue_map[[1]] %>% data.table()
biosample_key = methods_temp$biosample_map[[1]] %>% data.table()
key = paste0(GTEx_key[[1]], ".", biosample_key[[1]])
temp$pred_name_long = methods_temp$pred_name_long[1]
temp = filter(temp, tissue.biosample %in% key)
enr.all = temp

# loop through rest, filtering to tissue/biosample matches
if (length(methods) > 1) {
  for (i in 1:length(methods)) {
    for (j in 1:length(distances)){
    temp = read.table(file = enrTableFiles[1+(i-1)*length(distances)+j-1],header = TRUE, fill = TRUE) %>% drop_na()
    temp$method = methods[i]
    temp$distance = distances[j]
    temp$tissue.biosample = paste0(temp$GTExTissue, ".", temp$Biosample)
    methods_temp = dplyr::filter(methods_config, method==methods[i])
    GTEx_key = methods_temp$GTExTissue_map[[1]] %>% data.table()
    biosample_key = methods_temp$biosample_map[[1]] %>% data.table()
    key = paste0(GTEx_key[[1]], ".", biosample_key[[1]])
    temp$pred_name_long = methods_temp$pred_name_long[1]
    temp = filter(temp, tissue.biosample %in% key)
    enr.all = rbind(enr.all, temp)
    }
  }
}
enr.all = distinct(enr.all)

# get stats, edit names, define labels 
enr.all = enr.all[order(enr.all$enrichment), ]
max.enr = max(enr.all$enrichment)
enrLabel = 'Enrichment\n(GTEx variants/all common variants)'

enr.all$method = factor(enr.all$method, levels=methods, ordered=TRUE)

## aggregate and process prediction metrics
# predTable = os.path.join(config["outDir"], "{method}", "predictionTables", "GTExTissue{GTExTissue}.Biosample{Biosample}.byDistance.tsv")
# first one
method.this = methods[1]
methods_temp = dplyr::filter(methods_config, method==method.this)
GTEx_key = methods_temp$GTExTissue_map[[1]] %>% data.table()
biosample_key = methods_temp$biosample_map[[1]] %>% data.table()
key_file = paste0("GTExTissue", GTEx_key[[1]], ".", "Biosample", biosample_key[[1]], ".byDistance.tsv")
key_path = file.path(outDir, method.this, key_file)
temp = read.table(file = key_paths, header = TRUE, fill = TRUE) %>% drop_na()
temp$GTExTissue = GTEx_key[[1]]
temp$Biosample = biosample_key[[1]]
temp$method = methods.this
pred.all = temp

for (i in 1:length(methods)){
  method.this = methods[i]
  methods_temp = dplyr::filter(methods_config, method==method.this)
  GTEx_key = methods_temp$GTExTissue_map[[i]] %>% data.table()
  biosample_key = methods_temp$biosample_map[[1]] %>% data.table()
    for (j in 1:length(GTEx_key)){
      key_file = paste0("GTExTissue", GTEx_key[[j]], ".", "Biosample", biosample_key[[j]], ".byDistance.tsv")
      key_path = file.path(outDir, method.this, key_files)
      temp = read.table(file = key_path, header = TRUE, fill = TRUE) %>% drop_na()
      temp$GTExTissue = GTEx_key[[j]]
      temp$Biosample = biosample_key[[j]]
      temp$method = methods.this
      pred.all = rbind(pred.all, temp)
    }
}
pred.all = distinct(pred.all)

# add column with plotting name by left_joining with methods_config selecting just the method and plotting name columns
methods_config.select = dplyr::select(methods_config, method, pred_name_long)
pred.all = left_join(pred.all, methods_config.select, by=c("method"))
enr.all = left_join(enr.all, methods_config.select, by=c("method"))

# make y-axis labels (distance range, min-max variants)
pred.all$distance.label = 0
enr.all$distance.label = 0

distances_plus0 = c(0, distances_plus)
for (i in 2:length(distances_plus0)){
  distance.this = distances_plus0[i]
  if (distance.this == 30000000){
    label = "All variants\n"
  } else {
    label = paste0(distance_plus0[i-1], " - ", distance.this, "bp\n")
  }
  pred.filt = dplyr::filter(pred.all, distance==distance.this)
  label.full = paste0(label, "N = ", min(pred.filt$total.variants), " - ", max(pred.filt$total.variants), " variants")
  pred.all$distance.label[pred.all$distance==distance.this] = label.full
  enr.all$distance.label[enr.all$distance==distance.this] = label.full
}

write.table(pred.all, outPredMetrics, col.names=TRUE, quote=FALSE, sep="\t")
write.table(enr.all, outEnrAll, col.names=TRUE, quote=FALSE, sep="\t")


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
sr.predicted = ggplot(pred.all, aes(x = distance.label, y = correctGene.ifOverlaps, fill=pred_name_long)) +
  geom_boxplot() +
  scale_fill_manual(values=cpList) +
  theme_minimal() + ggtitle('Variants linked to correct gene,\ngiven overlapping predicted enhancer') + xlab('') +
  ylab('Fraction of GTEx variants\noverlapping predicted enhancers') +
  theme(axis.text.y = element_blank()) + coord_flip()

## save final plots
all.tissues =  ggarrange(enr.boxplot, sr.overlaps, sr.predicted, nrow=1, ncol=3)
pdf(file=outFile, width=12, height=5)
  print(all.tissues)
dev.off()
 
 

