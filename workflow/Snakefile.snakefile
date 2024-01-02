# coding: utf-8

import pandas as pd
from os.path import join

conda: "mamba"

# load prediction method config file
methods_config_file = config["methodsTable"]
methods_config = pd.read_table(methods_config_file, na_values="").fillna("None").set_index("method", drop=False)

# import workflows
pd.set_option('display.max_columns', None)

include: "./scripts/add_biosamples_and_files_to_config.py"
include: "./scripts/process_biosample_tissue_maps.py"

# find set of all GTEx tissues with biosample matches across methods 
GTExTissues_all = []
for x in config["methods"]:
	GTExTissues_all.extend(methods_config.loc[x, "GTExTissue_map"])
	GTExTissues_matched = [*set(GTExTissues_all)]

include: "./rules/preprocessing.smk"
include: "./rules/recall.smk"
include: "./rules/enrichment.smk"
include: "./rules/visualization.smk"

# generate files that need looping
variantPredictionsIntFiles = []
commonVarPredictionsIntFiles = []
enrichmentMatrices_threshold_Files = []
enrichmentMatrices_distance_Files = []
predTables_distance_Files = []
predTables_threshold_Files = []
enrichmentRecallTableFiles = []


for x in config["methods"]:
	variantPredictionsIntFiles.extend(expand(os.path.join(config["outDir"], x, "{biosample}", "GTExVariants-enhancerPredictionsInt.tsv.gz"), biosample=methods_config.loc[x, 'biosamples']))
	commonVarPredictionsIntFiles.extend(expand(os.path.join(config["outDir"], x, "{biosample}", "distalNoncodingBackgroundSNPs-enhancerPredictionsInt.tsv.gz"), biosample=methods_config.loc[x, 'biosamples']))
	enrichmentMatrices_distance_Files.extend(expand(os.path.join(config["outDir"], x, "enrichmentTables", "enrichmentTable.under{distance}bp.tsv"), distance=config["distances"]))
	predTables_distance_Files.extend(expand(os.path.join(config["outDir"], x, "predictionTables", "GTExTissue{GTExTissue}.Biosample{Biosample}.byDistance.tsv"), zip, GTExTissue=methods_config.loc[x, "GTExTissue_map"], Biosample=methods_config.loc[x, "biosample_map"]))
	predTables_threshold_Files.extend(expand(os.path.join(config["outDir"], x, "predictionTables", "GTExTissue{GTExTissue}.Biosample{Biosample}.byThreshold.tsv"), zip, GTExTissue=methods_config.loc[x, "GTExTissue_map"], Biosample=methods_config.loc[x, "biosample_map"]))
	enrichmentRecallTableFiles.extend(expand(os.path.join(config["outDir"], x, "enrichmentRecallTables", "GTExTissue{GTExTissue}.Biosample{Biosample}.tsv"), zip, GTExTissue=methods_config.loc[x, "GTExTissue_map"], Biosample=methods_config.loc[x, "biosample_map"]))

rule first:
	input:
		variantsPredictionsInt = variantPredictionsIntFiles,
		commonVarPredictionsInt = commonVarPredictionsIntFiles,
		variantsPerGTExTissue = expand(os.path.join(config["outDir"], "{method}", "nVariantsPerGTExTissue.tsv"), method=config["methods"])

rule second:
	input:
		thresholdSpans = expand(os.path.join(config["outDir"], "{method}", "thresholdSpan.tsv"), method=config["methods"]),
		enrichmentMatrices_distance = enrichmentMatrices_distance_Files,
		enrichmentTable = expand(os.path.join(config["outDir"], "{method}", "enrichmentTables", "giant_enrichmentTable.threshold.tsv"),method=config["methods"]),
		predTables_distance = predTables_distance_Files,
		predTables_threshold = predTables_threshold_Files

rule third:
	input:
		enrichmentRecallTables = enrichmentRecallTableFiles

rule fourth:
	input:
		enrichmentRecallCurves = expand(os.path.join(config["outDir"],  "plots", "enrichmentRecall.GTExTissue{GTExTissue}.pdf"), GTExTissue=GTExTissues_matched),
		aggregateEnrichmentRecall = os.path.join(config["outDir"], "plots", "aggregate_enrichment_recall.pdf"),
		enrichmentWithCI = expand(os.path.join(config["outDir"],  "plots", "enrichments.Recall{recall}.GTExTissue{GTExTissue}.tsv"), recall=config["recalls"], GTExTissue=GTExTissues_matched),
		comparisonPlots =  os.path.join(config["outDir"], "plots", "final_comparison_figure_grouped_by_dist.pdf"),

rule fifth:
	input:
		benchmarkingReport = os.path.join(config["outDir"], "benchmarking_report.html")

