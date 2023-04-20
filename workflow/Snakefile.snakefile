# coding: utf-8

import pandas as pd
from os.path import join

# config file containing samples & parameters
# configfile: "../config/config-hg38-v2.yml"

# load prediction method config file
methods_config_file = config["methodsTable"]
methods_config = pd.read_table(methods_config_file, na_values="").fillna("None").set_index("method", drop=False)
#methods_config["GTExTissue_map"] = methods_config["GTExTissue_map"].apply(eval)
#methods_config["biosample_map"] = methods_config["biosample_map"].apply(eval)

# import workflows
include: "./rules/preprocessing.smk"
include: "./rules/eQTL_predictions.smk"
include: "./rules/enrichment.smk"
include: "./rules/eQTL_thresholding.smk"
include: "./scripts/add_biosamples_and_files_to_config.py"
include:  "./scripts/generate_threshold_span.py"
include: "./scripts/process_biosample_tissue_maps.py"

# generate files that need looping
variantPredictionsIntFiles = []
commonVarPredictionsIntFiles = []
enrichmentMatrices_threshold_Files = []
predTables_distance_Files = []
predTables_threshold_Files = []

for x in config["methods"]:
	variantPredictionsIntFiles.extend(expand(os.path.join(config["outDir"], x, "{biosample}", "GTExVariants-enhancerPredictionsInt.tsv.gz"), biosample=methods_config.loc[x, 'biosamples']))
	commonVarPredictionsIntFiles.extend(expand(os.path.join(config["outDir"], x, "{biosample}", "distalNoncodingBackgroundSNPs-enhancerPredictionsInt.tsv.gz"), biosample=methods_config.loc[x, 'biosamples']))
	enrichmentMatrices_threshold_Files.extend(expand(os.path.join(config["outDir"], "{method}", "enrichmentTable.thresh{threshold}.tsv"), threshold=methods_config.loc[x, "thresholdSpan"]))
	predTables_distance_Files.extend(expand(os.path.join(config["outDir"], "{method}", predictionTables, "GTExTissue{GTExTissue}.Biosample{Biosample}.byDistance.tsv"), zip, GTExTissue=methods_config.loc[x, "GTExTissue_map"], Biosample=methods_config.loc[x, "biosample_map"]))
	predTables_threshold_Files.extend(expand(os.path.join(config["outDir"], "{method}", predictionTables, "GTExTissue{GTExTissue}.Biosample{Biosample}.byThreshold.tsv"), zip, GTExTissue=methods_config.loc[x, "GTExTissue_map"], Biosample=methods_config.loc[x, "biosample_map"]))

rule first:
	input:
		variantsPredictionsInt = variantPredictionsIntFiles,
		commonVarPredictionsInt = commonVarPredictionsIntFiles

rule second:
	input:
		enrichmentMatrices_distance = expand(os.path.join(config["outDir"], "{method}", "enrichmentTable.under{distance}bp.tsv"), method=config['methods'], distance=config["distances"]),
		enrichmentMatrices_threshold = enrichmentMatrices_thresholdFiles,
		predTables_distance = predTables_distance_Files,
		predTables_threshold = predTables_threshold_Files






########################### OLD ###############################

# generate output files
variantsByTissueFiles = []
predictionsByBiosampleFiles = []
predTablesFiles = []
thresholdTableFiles = []
varIntFiles = []
commonVarFiles = []
predThresholdedFiles = []
enrichmentTableFiles = []
countMatrixFiles = []

for x in config["methods"]:
	variantsByTissueFiles.extend(expand(os.path.join(config["outDir"], x, "eGenePrediction",  "{GTExTissue}.{biosample}.filteredVariants.tsv"), zip, GTExTissue=methods_config.loc[x, "GTExTissue_map"], biosample=methods_config.loc[x, "biosample_map"]))
	predictionsByBiosampleFiles.extend(expand(os.path.join(config["outDir"], x, "{biosample}", "enhancerPredictions.sorted.bed.gz"), biosample=methods_config.loc[x, "biosample_map"]))
	predTablesFiles.extend(expand(os.path.join(config["outDir"], x, "eGenePrediction", "{GTExTissue}.{Biosample}.predictionTable.tsv"), zip, GTExTissue=methods_config.loc[x, "GTExTissue_map"], Biosample=methods_config.loc[x, "biosample_map"]))
	varIntFiles.extend(expand(os.path.join(config["outDir"], x, "{biosample}", "GTExVariants-enhancerPredictionsInt.{threshold}.tsv.gz"), biosample=methods_config.loc[x, "biosamples"], threshold=methods_config.loc[x, "thresholdSpan"]))		
	commonVarFiles.extend(expand(os.path.join(config["outDir"], x, "{biosample}", "commonVarPerBiosample.{threshold}.tsv"), biosample=methods_config.loc[x, "biosamples"], threshold=methods_config.loc[x, "thresholdSpan"]))	
	predThresholdedFiles.extend(expand(os.path.join(config["outDir"], x, "{biosample}", "enhancerPredictions.thresholded.{threshold}.bed.gz"), biosample=methods_config.loc[x, "biosamples"], threshold=methods_config.loc[x, "thresholdSpan"]))
	thresholdTableFiles.extend(expand(os.path.join(config["outDir"], "thresholdTables", x, "{GTExTissue}.{biosample}.tsv"), zip, GTExTissue=methods_config.loc[x, "GTExTissue_map"], biosample=methods_config.loc[x, "biosample_map"]))
	enrichmentTableFiles.extend(expand(os.path.join(config["outDir"], "{method}", "enrichmentTable.{threshold}.tsv"), method=config["methods"], threshold=methods_config.loc[x, "thresholdSpan"]))
	countMatrixFiles.extend(expand(os.path.join(config["outDir"], "{method}", "count_matrix.{threshold}.tsv"), method=config["methods"], threshold=methods_config.loc[x, "thresholdSpan"]))

rule first:
	input:
		geneUniverses = expand(os.path.join(config["outDir"], "{method}", "geneUniverse.tsv"), method=config["methods"]),
		GTExExpressedGenes = os.path.join(config["outDir"], "generalReference", "GTExVariants.filtered.PIP0.5.distalNoncoding.expressed.tsv"),
		predictionsThresholded = predThresholdedFiles,
		predTablesForCalcs = predTablesFiles,
		samples = expand(os.path.join(config["outDir"], "{method}", "biosampleList.tsv"), method=config["methods"]),
		filteredGTExVar = os.path.join(config["outDir"], "generalReference/GTExVariants.filtered.PIP0.5.distalNoncoding.tsv.gz"),
		partitionDistalNoncoding = os.path.join(config["outDir"], "generalReference", "Partition.distalNoncoding.bed"),
		commonVarDistalNoncoding = os.path.join(config["outDir"], "generalReference", "distalNoncoding.bg.SNPs.bed.gz"),
		filteredGTExVariantsFinal = expand(os.path.join(config["outDir"], "{method}", "GTExVariants.filteredForMethod.tsv"), method=config["methods"]),
		variantsPredictionsInt = varIntFiles,
		variantsPerGTExTissue = expand(os.path.join(config["outDir"], "{method}", "nVariantsPerGTExTissue.tsv"), method=config["methods"]),
		commonVarPerBiosample = commonVarFiles,
		
rule second: 
	input:
		countMatrix = countMatrixFiles
		
rule third: 
	input:
		enrichmentTable = enrichmentTableFiles

rule fourth:
	input:
		variantsByTissue = variantsByTissueFiles,
		predictionsByBiosample = predictionsByBiosampleFiles,
		predTables = predTablesFiles,
		thresholdTables = thresholdTableFiles
		
rule fifth:
	input:
		ERCurveMeanFull = expand(os.path.join(config["outDir"], "{method}", "ERCurveMean.full.pdf"), method=config["methods"] ),
		ERCurveMeanZoom = expand(os.path.join(config["outDir"], "{method}", "ERCurveMean.zoom.pdf"), method=config["methods"]),
		ERCurveMaxFull = expand(os.path.join(config["outDir"], "{method}", "ERCurveMax.full.pdf"), method=config["methods"]),
		ERCurveMaxZoom = expand(os.path.join(config["outDir"], "{method}", "ERCurveMax.zoom.pdf"), method=config["methods"]),
		ERCurveTable = expand(os.path.join(config["outDir"], "{method}", "ERCurveTable.tsv"), method=config["methods"])
		
################################################################################################################################

