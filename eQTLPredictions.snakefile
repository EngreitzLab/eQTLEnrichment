# coding: utf-8

import pandas as pd
from os.path import join

# import workflows
include: "./rules/eQTL_predictions.smk"
include: "./rules/eQTL_enrichment.smk"

# config file containing samples & parameters
configfile: "./config/config.yml"

# generate files across methods
filteredPredictions = []
predTables = []

for x in config["methods"]:

	# list of prediction tables
	predTables.extend(expand(os.path.join(config["outDir"], x, "{GTExTissue}.{Biosample}.predictionTable.tsv"), zip, GTExTissue=config["predRates"]["mapGTExTissues"][x], Biosample=config["predRates"]["mapBiosamples"][x]))
	
	# list of filtered enhancer-gene predictions
	filteredPredictions.extend(expand(os.path.join(config["outDir"], x, "{Biosample}.filteredPredictions.tsv.gz"), Biosample=config["predRates"]["mapBiosamples"][x]))
	
	# list of output sensitivity plot names
	sensitivityPlots.extend(expand(os.path.join(config["outDir"], "{GTExTissue}.{method}Enhancers.variantOverlapSensitivity.pdf"), GTExTissue=config["predRates"]["mapGTExTissues"][x], method=x))

rule all:
	input:
		sortedPredictions = expand(os.path.join(config["outDir"], "{method}/enhancerPredictions.sorted.bed.gz"), method=config["methods"]),
		samples = expand(os.path.join(config["outDir"], "{method}/biosampleList.tsv"), method=config["methods"]),
		filteredVariantsPredictionsInt = expand(os.path.join(config["outDir"], "{method}/GTExVariants-enhancerPredictions.PIP0.5.distalNoncoding.tsv.gz"), method=config["methods"]),
		commonVariantsPredictionsInt = expand(os.path.join(config["outDir"], "{method}/commonVariants-enhancerPredictions.distalNoncoding.tsv.gz"), method=config["methods"]),
		filteredVariantsEnrichment = os.path.join(config["outDir"], "variantFilesForEnrichment/GTEx.filtered.PIP0.5.tsv.gz"),
		partitionDistalNoncoding = os.path.join(config["outDir"], "variantFilesForEnrichment/Partition.distalNoncoding.bed"),
		filteredVariantsEnrichmentDistalNoncoding = os.path.join(config["outDir"], "variantFilesForEnrichment/GTEx.filtered.PIP0.5.distalNoncoding.tsv.gz"),
		commonVariantsDistalNoncoding = os.path.join(config["outDir"], "variantFilesForEnrichment/distalNoncoding.bg.SNPs.bed.gz"),
		variantsPerGTExTissue = os.path.join(config["outDir"], "variantFilesForEnrichment/variantsPerGTExTissue.tsv"),
		enrichmentTables = expand(os.path.join(config["outDir"], "{method}/enrichmentTable.tsv"), method=config["methods"]),
		countMatrices = expand(os.path.join(config["outDir"], "{method}/countMatrix.tsv"), method=config["methods"]),
		enhancerSetSizes = expand(os.path.join(config["outDir"], "{method}/basesPerEnhancerSet.tsv"), method=config["methods"]),
		commonVariantsPerBiosample = expand(os.path.join(config["outDir"], "{method}/commonVarPerBiosample.tsv"), method=config["methods"]),
		fullEnrichmentHeatmaps = expand(os.path.join(config["outDir"], "{method}/enrichmentHeatmap.full.pdf"), method=config["methods"]),
		aggregatedEnrichmentHeatmaps = expand(os.path.join(config["outDir"], "{method}/enrichmentHeatmap.aggregated.pdf"), method=config["methods"]),
		cdfPlot = os.path.join(config["outDir"], "cdf.pdf"),
		densityPlot = os.path.join(config["outDir"], "density.pdf"),
		filteredVariantsPrediction = expand(os.path.join(config["outDir"], "variantFilesForPrediction", "{GTExTissue}.filteredVariants.tsv.gz"), GTExTissue=config["GTExTissues"]),
		filteredVariantsProximal = expand(os.path.join(config["outDir"], "variantFilesForPrediction", "{GTExTissue}.filteredVariants.proximalGenes.tsv.gz"), GTExTissue=config["GTExTissues"]),
		allFilteredPredictions = filteredPredictions,
		allPredictionTables = predTables,
		allSensitivityPlots = sensitivityPlots,
		predictionPlot = os.path.join(config["outDir"], "predictionRates.pdf")