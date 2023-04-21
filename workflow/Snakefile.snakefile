# coding: utf-8

import pandas as pd
from os.path import join

# config file containing samples & parameters
# configfile: "../config/config-hg38-v2.yml"

# load prediction method config file
methods_config_file = config["methodsTable"]
methods_config = pd.read_table(methods_config_file, na_values="").fillna("None").set_index("method", drop=False)


# import workflows
pd.set_option('display.max_columns', None)

include: "./scripts/add_biosamples_and_files_to_config.py"
include: "./scripts/generate_threshold_span.py"
include: "./scripts/process_biosample_tissue_maps.py"

include: "./rules/preprocessing.smk"
include: "./rules/recall.smk"
include: "./rules/enrichment.smk"

# generate files that need looping
variantPredictionsIntFiles = []
commonVarPredictionsIntFiles = []
enrichmentMatrices_threshold_Files = []
predTables_distance_Files = []
predTables_threshold_Files = []

for x in config["methods"]:
	variantPredictionsIntFiles.extend(expand(os.path.join(config["outDir"], x, "{biosample}", "GTExVariants-enhancerPredictionsInt.tsv.gz"), biosample=methods_config.loc[x, 'biosamples']))
	commonVarPredictionsIntFiles.extend(expand(os.path.join(config["outDir"], x, "{biosample}", "distalNoncodingBackgroundSNPs-enhancerPredictionsInt.tsv.gz"), biosample=methods_config.loc[x, 'biosamples']))
	enrichmentMatrices_threshold_Files.extend(expand(os.path.join(config["outDir"], x, "enrichmentTable.thresh{threshold}.tsv"), threshold=methods_config.loc[x, "thresholdSpan"]))
	predTables_distance_Files.extend(expand(os.path.join(config["outDir"], x, "predictionTables", "GTExTissue{GTExTissue}.Biosample{Biosample}.byDistance.tsv"), zip, GTExTissue=methods_config.loc[x, "GTExTissue_map"], Biosample=methods_config.loc[x, "biosample_map"]))
	predTables_threshold_Files.extend(expand(os.path.join(config["outDir"], x, "predictionTables", "GTExTissue{GTExTissue}.Biosample{Biosample}.byThreshold.tsv"), zip, GTExTissue=methods_config.loc[x, "GTExTissue_map"], Biosample=methods_config.loc[x, "biosample_map"]))

rule first:
	input:
		variantsPredictionsInt = variantPredictionsIntFiles,
		commonVarPredictionsInt = commonVarPredictionsIntFiles,
		variantsPerGTExTissue = expand(os.path.join(config["outDir"], "{method}", "nVariantsPerGTExTissue.tsv"), method=config["methods"])

rule second:
	input:
		enrichmentMatrices_distance = expand(os.path.join(config["outDir"], "{method}", "enrichmentTable.under{distance}bp.tsv"), method=config['methods'], distance=config["distances"]),
		enrichmentMatrices_threshold = enrichmentMatrices_threshold_Files,
		predTables_distance = predTables_distance_Files,
		predTables_threshold = predTables_threshold_Files
