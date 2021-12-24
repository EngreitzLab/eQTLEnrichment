# coding: utf-8

import pandas as pd
from os.path import join

# config file containing samples & parameters
# configfile: "../config/config-hg38-v2.yml"

# load prediction method config file
methods_config_file = config["methodsTable"]
methods_config = pd.read_table(methods_config_file, na_values="").fillna("None").set_index("method", drop=False)
methods_config["GTExTissue_map"] = methods_config["GTExTissue_map"].apply(eval)
methods_config["biosample_map"] = methods_config["biosample_map"].apply(eval)

# import workflows
include: "./rules/eQTL_predictions.smk"
include: "./rules/eQTL_enrichment.smk"
include: "./scripts/add_biosamples_and_files_to_config.py"

# generate output files
variantsByTissueFiles = []
variantsByTissueProximalFiles = []
predictionsByBiosampleFiles = []
predTablesForCalcsFiles = []
predTablesFiles = []
thresholdTableFiles = []
sensitivityPlotsFiles = []
varIntFiles = []
commonVarFiles = []
predThresholdedFiles = []
basesPerSetFiles = []

for x in config["methods"]:
	# list of variant files
	variantsByTissueFiles.extend(expand(os.path.join(config["outDir"], x, "eGenePrediction",  "{GTExTissue}.{biosample}.filteredVariants.tsv"), zip, GTExTissue=methods_config.loc[x, "GTExTissue_map"], biosample=methods_config.loc[x, "biosample_map"]))
	# list of variant files with proximal genes
	variantsByTissueProximalFiles.extend(expand(os.path.join(config["outDir"], x, "eGenePrediction",  "{GTExTissue}.{biosample}.filteredVariants.proximal.tsv"), zip, GTExTissue=methods_config.loc[x, "GTExTissue_map"], biosample=methods_config.loc[x, "biosample_map"]))
	# list of filtered prediction files
	predictionsByBiosampleFiles.extend(expand(os.path.join(config["outDir"], x, "{biosample}", "enhancerPredictions.sorted.bed.gz"), biosample=methods_config.loc[x, "biosample_map"]))
	# list of prediction tables
	predTablesFiles.extend(expand(os.path.join(config["outDir"], x, "eGenePrediction", "{GTExTissue}.{Biosample}.predictionTable.tsv"), zip, GTExTissue=methods_config.loc[x, "GTExTissue_map"], Biosample=methods_config.loc[x, "biosample_map"]))
	predTablesForCalcsFiles.extend(expand(os.path.join(config["outDir"], x, "eGenePrediction", "{GTExTissue}.{Biosample}.predictionTable.forCalcs.tsv"), zip, GTExTissue=methods_config.loc[x, "GTExTissue_map"], Biosample=methods_config.loc[x, "biosample_map"]))
	# list of output plot names
	sensitivityPlotsFiles.extend(expand(os.path.join(config["outDir"], "sensitivityPlots", "{GTExTissue}.{method}Enhancers.variantOverlapSensitivity.pdf"), GTExTissue=methods_config.loc[x, "GTExTissue_map"], method=x))
	# list of variant-prediction intersections
	varIntFiles.extend(expand(os.path.join(config["outDir"], x, "{biosample}", "GTExVariants-enhancerPredictionsInt.tsv.gz"), biosample=methods_config.loc[x, "biosamples"]))		
	commonVarFiles.extend(expand(os.path.join(config["outDir"], x, "{biosample}", "commonVarPerBiosample.tsv"), biosample=methods_config.loc[x, "biosamples"]))	
	predThresholdedFiles.extend(expand(os.path.join(config["outDir"], x, "{biosample}", "enhancerPredictions.thresholded.bed.gz"), biosample=methods_config.loc[x, "biosamples"]))
	basesPerSetFiles.extend(expand(os.path.join(config["outDir"], x, "biosample", "basesPerEnhancerSet.tsv"), biosample=methods_config.loc[x, "biosamples"]))	
	thresholdTableFiles.extend(expand(os.path.join(config["outDir"], "thresholdTables", x, "{GTExTissue}.{Biosample}.tsv"), zip, GTExTissue=methods_config.loc[x, "GTExTissue_map"], Biosample=methods_config.loc[x, "biosample_map"]))

rule first:
	input:
		geneUniverses = expand(os.path.join(config["outDir"], "{method}", "geneUniverse.tsv"), method=config["methods"]),
		GTExExpressedGenes = os.path.join(config["outDir"], "generalReference", "GTExVariants.filtered.PIP0.5.distalNoncoding.expressed.tsv"),
		predictionsThresholded = predThresholdedFiles,
		predTablesForCalcs = predTablesForCalcsFiles,
		samples = expand(os.path.join(config["outDir"], "{method}", "biosampleList.tsv"), method=config["methods"]),
		filteredGTExVar = os.path.join(config["outDir"], "generalReference/GTExVariants.filtered.PIP0.5.distalNoncoding.tsv.gz"),
		partitionDistalNoncoding = os.path.join(config["outDir"], "generalReference", "Partition.distalNoncoding.bed"),
		commonVarDistalNoncoding = os.path.join(config["outDir"], "generalReference", "distalNoncoding.bg.SNPs.bed.gz"),
		filteredGTExVariantsFinal = expand(os.path.join(config["outDir"], "{method}", "GTExVariants.filteredForMethod.tsv"), method=config["methods"]),
		variantsPredictionsInt = varIntFiles,
		variantsPerGTExTissue = expand(os.path.join(config["outDir"], "{method}", "nVariantsPerGTExTissue.tsv"), method=config["methods"]),
		commonVarPerBiosample = commonVarFiles,
		#basesPerEnhancerSet=basesPerSetFiles
		
rule second: 
	input:
		countMatrix = expand(os.path.join(config["outDir"], "{method}", "count_matrix.tsv"), method=config["methods"]),

rule third: 
	input:
		enrichmentTable = expand(os.path.join(config["outDir"], "{method}", "enrichmentTable.tsv"), method=config["methods"]),
		cdf = os.path.join(config["outDir"], "cdf.pdf"),
		density = os.path.join(config["outDir"], "density.pdf"),
		boxplot = os.path.join(config["outDir"], "boxplot.pdf"),
		colorPalette = os.path.join(config["outDir"], "colorPalette.rds")
		

rule fourth:
	input:
		variantsByTissue = variantsByTissueFiles,
		predictionsByBiosample = predictionsByBiosampleFiles,
		predTables = predTablesFiles,
		thresholdTables = thresholdTableFiles,
		predictionPlot = os.path.join(config["outDir"], "predictionRates.pdf"),
		PPVPlot = os.path.join(config["outDir"], "PPV.pdf"),
		predictionMetrics = os.path.join(config["outDir"],"predictionMetrics.tsv"),
		#sensitivitiesTable = os.path.join(config["outDir"],"sensitivitiesTable.tsv"),
		enrichmentReport = os.path.join(config["outDir"], "enrichment_report.html")
		
rule ER_only:
	input:
		cdf = os.path.join(config["outDir"], "cdf.pdf"),
		density = os.path.join(config["outDir"], "density.pdf"),
		boxplot = os.path.join(config["outDir"], "boxplot.pdf"),
		colorPalette = os.path.join(config["outDir"], "colorPalette.rds"),
		enrichmentReport = os.path.join(config["outDir"], "enrichment_report.html")
	
		
################################################################################################################################

