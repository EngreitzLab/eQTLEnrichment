
rule compute_prediction_table_by_distance:
	input:
		variantsPredictionsInt = os.path.join(config["outDir"], "{method}", "biosamples", "{Biosample}", "GTExVariants-enhancerPredictionsInt.tsv.gz"),
		filteredGTExVariantsFinal = os.path.join(config["outDir"], "{method}", "intermediate", "GTExVariants.filteredForMethod.tsv.gz")
	params:
		distances_max = config["distances_max"],
		distances_min = config["distances_min"],
		threshold = lambda wildcards: methods_config.loc[wildcards.method, "threshold"],
	output:
		predTable = os.path.join(config["outDir"], "{method}", "predictionTables", "GTExTissue{GTExTissue}.Biosample{Biosample}.byDistance.tsv")
	resources:
		mem_mb = determine_mem_mb
	conda:
		os.path.join(config["envDir"], "eQTLEnv.yml")
	script:
		os.path.join(config["codeDir"], "enrichment_recall", "compute_prediction_table_distance.R")


rule compute_prediction_table_by_threshold:
	input:
		variantsPredictionsInt = os.path.join(config["outDir"], "{method}", "biosamples", "{Biosample}", "GTExVariants-enhancerPredictionsInt.tsv.gz"),
		filteredGTExVariantsFinal = os.path.join(config["outDir"], "{method}", "intermediate", "GTExVariants.filteredForMethod.tsv.gz"),
		thresholdSpan = os.path.join(config["outDir"], "{method}", "intermediate", "thresholdSpan.tsv")
	output:
		predTable = os.path.join(config["outDir"], "{method}", "predictionTables", "GTExTissue{GTExTissue}.Biosample{Biosample}.byThreshold.tsv")
	resources:
		mem_mb = determine_mem_mb
	conda:
		os.path.join(config["envDir"], "eQTLEnv.yml")
	script:
		os.path.join(config["codeDir"], "enrichment_recall", "compute_prediction_table_threshold.R")


