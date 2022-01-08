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
include:  "./scripts/generate_threshold_span.py"

# generate output files
variantsByTissueFiles = []
predictionsByBiosampleFiles = []
predTablesForCalcsFiles = []
predTablesFiles = []
thresholdTableFiles = []
sensitivityPlotsFiles = []
varIntFiles = []
commonVarFiles = []
predThresholdedFiles = []
enrichmentTableFiles = []
countMatrixFiles = []

for x in config["methods"]:
	# list of variant files
	variantsByTissueFiles.extend(expand(os.path.join(config["outDir"], x, "eGenePrediction",  "{GTExTissue}.{biosample}.filteredVariants.tsv"), zip, GTExTissue=methods_config.loc[x, "GTExTissue_map"], biosample=methods_config.loc[x, "biosample_map"]))
	# list of filtered prediction files
	predictionsByBiosampleFiles.extend(expand(os.path.join(config["outDir"], x, "{biosample}", "enhancerPredictions.sorted.bed.gz"), biosample=methods_config.loc[x, "biosample_map"]))
	# list of prediction tables
	predTablesFiles.extend(expand(os.path.join(config["outDir"], x, "eGenePrediction", "{GTExTissue}.{Biosample}.predictionTable.tsv"), zip, GTExTissue=methods_config.loc[x, "GTExTissue_map"], Biosample=methods_config.loc[x, "biosample_map"]))
	predTablesForCalcsFiles.extend(expand(os.path.join(config["outDir"], x, "eGenePrediction", "{GTExTissue}.{Biosample}.predictionTable.forCalcs.tsv"), zip, GTExTissue=methods_config.loc[x, "GTExTissue_map"], Biosample=methods_config.loc[x, "biosample_map"]))
	# list of variant-prediction intersections
	varIntFiles.extend(expand(os.path.join(config["outDir"], x, "{biosample}", "GTExVariants-enhancerPredictionsInt.{threshold}.tsv.gz"), biosample=methods_config.loc[x, "biosamples"], threshold=methods_config.loc[x, "biosamples"]))		
	commonVarFiles.extend(expand(os.path.join(config["outDir"], x, "{biosample}", "commonVarPerBiosample.{threshold}.tsv"), biosample=methods_config.loc[x, "biosamples"], threshold=methods_config.loc[x, "thresholdSpan"]))	
	predThresholdedFiles.extend(expand(os.path.join(config["outDir"], x, "{biosample}", "enhancerPredictions.{threshold}.bed.gz"), biosample=methods_config.loc[x, "biosamples"], threshold=methods_config.loc[x, "thresholdSpan"]))
	thresholdTableFiles.extend(expand(os.path.join(config["outDir"], "thresholdTables", x, "{GTExTissue}.{biosample}.tsv"), zip, GTExTissue=methods_config.loc[x, "GTExTissue_map"], biosample=methods_config.loc[x, "biosample_map"]))
	enrichmentTablesFiles.extend(expand(os.path.join(config["outDir"], "{method}", "enrichmentTable.{threshold}.tsv"), method=config["methods"], threshold=methods_config.loc[x, "thresholdSpan"]))
	countMatrixFiles.extend(expand(os.path.join(config["outDir"], "{method}", "count_matrix.{threshold}.tsv"), method=config["methods"], threshold=methods_config.loc[x, "thresholdSpan"]))

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
		thresholdTables = thresholdTableFiles,
	
		
################################################################################################################################

