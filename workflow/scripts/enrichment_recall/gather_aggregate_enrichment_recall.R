suppressPackageStartupMessages({library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)})

# load data
predTableFiles = (snakemake@input$predTables)  %>% strsplit(" ") %>% unlist()
enrTableFile = (snakemake@input$enrichmentTable)
outERTable = (snakemake@output$ERCurveTable)
method.this = snakemake@wildcards$method

# read in prediction tables and aggregate
for (i in 1:length(predTableFiles)){
	predTable.this = read.table(predTableFiles[i], header=TRUE, sep="\t")
	predTable.this = dplyr::mutate(predTable.this, variants.linking.single = total.variants * recall.linking, variants.overlap.single = total.variants * recall.total, total.variants.single = total.variants) %>%
	dplyr::select(threshold, variants.linking.single, variants.overlap.single, total.variants.single)

	if (i==1){predTable.all = predTable.this} else {
		predTable.all = rbind(predTable.all, predTable.this)
	}	
}
predTable.all = group_by(predTable.all, threshold) %>%
	summarize(variants.linking = sum(variants.linking.single), variants.overlap = sum(variants.overlap.single), total.variants = sum(total.variants.single)) %>%
	mutate(recall.linking = variants.linking/total.variants, recall.total = variants.overlap/total.variants)

# read in enrichment table
enrTable = read.table(enrTableFile, header=TRUE, sep="\t")
this.enrTable = dplyr::filter(enrTable, GTExTissue=="all_matches", Biosample=="all_matches") %>% 
	dplyr::select(threshold, enrichment, CI_enr_low, CI_enr_high, SE_log_enr)

# merge by threshold with pred table
predTable.all = dplyr::left_join(predTable.all, this.enrTable, by="threshold")
predTable.all$method = method.this

# save table
write.table(predTable.all, outERTable, quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)

