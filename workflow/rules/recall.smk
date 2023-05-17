
rule compute_prediction_table_by_distance:
	input:
		variantsPredictionsInt = os.path.join(config["outDir"], "{method}", "{Biosample}", "GTExVariants-enhancerPredictionsInt.tsv.gz"),
		filteredGTExVariantsFinal = os.path.join(config["outDir"], "{method}", "GTExVariants.filteredForMethod.tsv")
	params:
		distances = config["distances"],
		threshold = lambda wildcards: methods_config.loc[wildcards.method, "threshold"],
	output:
		predTable = os.path.join(config["outDir"], "{method}", "predictionTables", "GTExTissue{GTExTissue}.Biosample{Biosample}.byDistance.tsv")
	conda:
			os.path.join(config["envDir"], "eQTLEnv.yml")
	script:
		os.path.join(config["codeDir"], "compute_prediction_table_distance.R")


rule compute_prediction_table_by_threshold:
	input:
		variantsPredictionsInt = os.path.join(config["outDir"], "{method}", "{Biosample}", "GTExVariants-enhancerPredictionsInt.tsv.gz"),
		filteredGTExVariantsFinal = os.path.join(config["outDir"], "{method}", "GTExVariants.filteredForMethod.tsv"),
		thresholdSpan = os.path.join(config["outDir"], "{method}", "thresholdSpan.tsv")
	params:
		thresholds = lambda wildcards: methods_config.loc[wildcards.method, "thresholdSpan"]
	output:
		predTable = os.path.join(config["outDir"], "{method}", "predictionTables", "GTExTissue{GTExTissue}.Biosample{Biosample}.byThreshold.tsv")
	conda:
			os.path.join(config["envDir"], "eQTLEnv.yml")
	script:
		os.path.join(config["codeDir"], "compute_prediction_table_threshold.R")


