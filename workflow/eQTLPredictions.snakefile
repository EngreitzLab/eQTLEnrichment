# coding: utf-8

import pandas as pd
from os.path import join

# import workflows
include: "./rules/eQTL_predictions.smk"
include: "./rules/eQTL_enrichment.smk"

# config file containing samples & parameters
# configfile: "../config/config.yml"

# generate output files
variantsByTissueFiles = []
variantsByTissueProximalFiles = []
predictionsByBiosampleFiles = []
predTablesFiles = []
sensitivityPlotsFiles = []

for x in config["methods"]:
	# list of variant files
	variantsByTissueFiles.extend(expand(os.path.join(config["outDir"], x, "eGenePrediction",  "{GTExTissue}.filteredVariants.tsv"), GTExTissue=config["predRates"]["mapGTExTissues"][x]))
	# list of variant files with proximal genes
	variantsByTissueProximalFiles.extend(expand(os.path.join(config["outDir"], x, "eGenePrediction",  "{GTExTissue}.filteredVariants.proximal.tsv"), GTExTissue=config["predRates"]["mapGTExTissues"][x]))
	# list of filtered prediction files
	predictionsByBiosampleFiles.extend(expand(os.path.join(config["outDir"], x, "eGenePrediction", "{Biosample}.filteredPredictions.tsv"), Biosample=config["predRates"]["mapBiosamples"][x]))
	# list of prediction tables
	predTablesFiles.extend(expand(os.path.join(config["outDir"], x, "eGenePrediction", "{GTExTissue}.{Biosample}.predictionTable.tsv"), zip, GTExTissue=config["predRates"]["mapGTExTissues"][x], Biosample=config["predRates"]["mapBiosamples"][x]))
	# list of output plot names
	sensitivityPlotsFiles.extend(expand(os.path.join(config["outDir"], "sensitivityPlots", "{GTExTissue}.{method}Enhancers.variantOverlapSensitivity.pdf"), GTExTissue=config["predRates"]["mapGTExTissues"][x], method=x))
	

rule all:
	input:
		geneUniverses = expand(os.path.join(config["outDir"], "{method}", "geneUniverse.tsv"), method=config["methods"]),
		GTExExpressedGenes = os.path.join(config["outDir"], "generalReference", "GTExVariants.filtered.PIP0.5.distalNoncoding.expressed.tsv"),
		predictionsSorted = expand(os.path.join(config["outDir"], "{method}", "enhancerPredictions.sorted.bed.gz"), method=config["methods"]),
		samples = expand(os.path.join(config["outDir"], "{method}", "biosampleList.tsv"), method=config["methods"]),
		filteredGTExVar = os.path.join(config["outDir"], "generalReference/GTExVariants.filtered.PIP0.5.distalNoncoding.tsv.gz"),
		partitionDistalNoncoding = os.path.join(config["outDir"], "generalReference", "Partition.distalNoncoding.bed"),
		commonVarDistalNoncoding = os.path.join(config["outDir"], "generalReference", "distalNoncoding.bg.SNPs.bed.gz"),
		filteredGTExVariantsFinal = expand(os.path.join(config["outDir"], "{method}", "GTExVariants.filteredForMethod.tsv"), method=config["methods"]),
		variantsPredictionsInt = expand(os.path.join(config["outDir"], "{method}", "GTExVariants-enhancerPredictionsInt.tsv.gz"), method=config["methods"]),
		countMatrix = expand(os.path.join(config["outDir"], "{method}", "countMatrix.tsv"), method=config["methods"]),
		variantsPerGTExTissue = expand(os.path.join(config["outDir"], "{method}", "nVariantsPerGTExTissue.tsv"), method=config["methods"]),
		commonVarPerBiosample = expand(os.path.join(config["outDir"], "{method}", "commonVarPerBiosample.tsv"), method=config["methods"]),
		enrichmentTable = expand(os.path.join(config["outDir"], "{method}", "enrichmentTable.tsv"), method=config["methods"]),
		#basesPerEnhancerSet = expand(os.path.join(config["outDir"], "{method}/basesPerEnhancerSet.tsv"), method=config["methods"]),
		#heatmapFull = expand(os.path.join(config["outDir"], "{method}", "enrichmentHeatmap.full.pdf"), method=config["methods"]),
		#heatmapAggregated = expand(os.path.join(config["outDir"], "{method}", "enrichmentHeatmap.aggregated.pdf"), method=config["methods"]),
		cdf = os.path.join(config["outDir"], "cdf.pdf"),
		density = os.path.join(config["outDir"], "density.pdf"),
		boxplot = os.path.join(config["outDir"], "boxplot.pdf"),
		enrichmentReport = os.path.join(config["outDir"], "enrichment_report.html"),
		variantsByTissue = variantsByTissueFiles,
		variantsByTissueProximal = variantsByTissueProximalFiles,
		predictionsByBiosample = predictionsByBiosampleFiles,
		predTables = predTablesFiles,
		sensitivityPlots = sensitivityPlotsFiles,
		predictionPlot = os.path.join(config["outDir"], "predictionRates.pdf"),
		PPVPlot = os.path.join(config["outDir"], "PPV.pdf"),
		predictionMetrics = os.path.join(config["outDir"],"predictionMetrics.tsv"),
		sensitivitiesTable = os.path.join(config["outDir"],"sensitivitiesTable.tsv")
		
################################################################################################################################

