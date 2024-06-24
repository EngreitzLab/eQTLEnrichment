
rule compute_prediction_table_by_distance:
	input:
		varPredInt = lambda wildcards: [os.path.join(config["outDir"], wildcards.method, "biosamples", biosample, "GTExVariants-enhancerPredictionsInt.tsv.gz") for biosample in methods_config.loc[wildcards.method, "biosamples"]],
		filteredGTExVariantsFinal = os.path.join(config["outDir"], "{method}", "intermediate", "GTExVariants.filteredForMethod.tsv.gz")
	params:
		distances_max = config["distances_max"],
		distances_min = config["distances_min"],
		threshold = lambda wildcards: methods_config.loc[wildcards.method, "threshold"],
		biosamples = lambda wildcards: methods_config.loc[wildcards.method, "biosamples"]
	output:
		predTable = os.path.join(config["outDir"], "{method}", "recallTables", "recallTable.byDistance.tsv")
	resources:
		mem_mb = determine_mem_mb
	conda:
		os.path.join(config["envDir"], "eQTLEnv.yml")
	script:
		os.path.join(config["codeDir"], "enrichment_recall", "compute_prediction_table_distance.R")


rule compute_prediction_table_by_threshold:
	input:
		varPredInt = lambda wildcards: [os.path.join(config["outDir"], wildcards.method, "biosamples", biosample, "GTExVariants-enhancerPredictionsInt.tsv.gz") for biosample in methods_config.loc[wildcards.method, "biosamples"]],
		filteredGTExVariantsFinal = os.path.join(config["outDir"], "{method}", "intermediate", "GTExVariants.filteredForMethod.tsv.gz"),
		thresholdSpan = os.path.join(config["outDir"], "{method}", "intermediate", "thresholdSpan.tsv"),
		map = os.path.join(config["outDir"], "{method}", "intermediate", "GTExTissueBiosampleMap.tsv")
	params:
		biosamples = lambda wildcards: methods_config.loc[wildcards.method, "biosamples"]
	output:
		predTable = os.path.join(config["outDir"], "{method}", "recallTables", "recallTable.acrossThresholds.tsv.gz")
	resources:
		mem_mb = determine_mem_mb
	conda:
		os.path.join(config["envDir"], "eQTLEnv.yml")
	script:
		os.path.join(config["codeDir"], "enrichment_recall", "compute_prediction_table_threshold.R")


