# color palette 	
rule generate_color_palette:
	params:
		user_inputs = [methods_config.loc[method, "color"] for method in config["methods"]],
		names = config["methods"],
		methods_config = config["methodsTable"]
	output:
		colorPalette = os.path.join(config["outDir"], "plots", "colorPalette.tsv"),		
	conda: 
		os.path.join(config["envDir"], "eQTLEnv.yml")
	script: 
		os.path.join(config["codeDir"], "visualization", "color_palette.R")

# gather data for enrichment recall curve per GTEx tissue/biosample match per method
rule gather_enrichment_recall:
	input: 
		predTable = os.path.join(config["outDir"], "{method}", "recallTables", "recallTable.acrossThresholds.tsv.gz"),
		enrichmentTable = os.path.join(config["outDir"], "{method}", "enrichmentTables", "enrichmentTable.acrossThresholds.tsv.gz"),
	conda: 
		os.path.join(config["envDir"], "eQTLEnv.yml")
	output:
		ERCurveTable = os.path.join(config["outDir"], "{method}", "enrichmentRecallTables", "GTExTissue{GTExTissue}.Biosample{biosample}.tsv")
	resources:
		mem_mb = determine_mem_mb
	script:
		os.path.join(config["codeDir"], "enrichment_recall", "gather_enrichment_recall.R")


# make table for "aggregate" enrichment/recall for all tissue/biosample matches per method
rule gather_enrichment_recall_aggregate:
	input: 
		predTable = os.path.join(config["outDir"], "{method}", "recallTables", "recallTable.acrossThresholds.tsv.gz"),
		enrichmentTable = os.path.join(config["outDir"], "{method}", "enrichmentTables", "enrichmentTable.acrossThresholds.tsv.gz"),
	conda: 
		os.path.join(config["envDir"], "eQTLEnv.yml")
	output:
		ERCurveTable = os.path.join(config["outDir"], "{method}", "enrichmentRecallTables", "AllMatches.tsv")
	resources:
		mem_mb = determine_mem_mb
	script:
		os.path.join(config["codeDir"], "enrichment_recall", "gather_aggregate_enrichment_recall.R")	

# combined enrichment recall curves-- one per GTEx tissue-- across methods
def get_table_files(GTExTissue):
	if (GTExTissue=="AllMatches"):
		tissue_files  = [os.path.join(config["outDir"], method, "enrichmentRecallTables", "AllMatches.tsv") for method in config["methods"]]
	else:
		substr = "GTExTissue" + GTExTissue
		all_files = flatten([[os.path.join(config["outDir"], method, "enrichmentRecallTables", f"GTExTissue{tissue}.Biosample{biosample}.tsv") for tissue, biosample in zip(methods_config.loc[method, "GTExTissue_map"], methods_config.loc[method, "biosample_map"])] for method in config['methods']])
		tissue_files =  [x for x in all_files if substr in x]
	return tissue_files

rule plot_enrichment_recall_curve:
	input:
		colorPalette = os.path.join(config["outDir"], "plots", "colorPalette.tsv"),
		enrichmentRecall_files = lambda wildcards: get_table_files(wildcards.GTExTissue)
	output:
		er_combined = os.path.join(config["outDir"],  "plots", "enrichmentRecall", "enrichmentRecall.GTExTissue{GTExTissue}.pdf"),
		er_combined_table = os.path.join(config["outDir"],  "plots", "enrichmentRecall", "enrichmentRecall.GTExTissue{GTExTissue}.tsv")
	resources:
		mem_mb = determine_mem_mb
	conda:
		os.path.join(config["envDir"], "eQTLEnv.yml")
	script:
		os.path.join(config["codeDir"], "visualization", "plot_enrichment_recall.R")
		
# plot enrichment with confidence intervals across methods, once per GTEx tissue and user-defined recalls
rule plot_enrichment_with_ci:
	input:
		colorPalette = os.path.join(config["outDir"], "plots", "colorPalette.tsv"),
		enrichmentRecall_files = lambda wildcards: get_table_files(wildcards.GTExTissue)
	params:
		thresholdPval = config["thresholdPval"]
	output:
		enr_at_recall = os.path.join(config["outDir"],  "plots", "enrichmentAtRecall", "enrichments.Recall{recall}.GTExTissue{GTExTissue}.pdf"),
		enr_at_recall_table = os.path.join(config["outDir"],  "plots", "enrichmentAtRecall", "enrichments.Recall{recall}.GTExTissue{GTExTissue}.tsv"),
		sign_table = os.path.join(config["outDir"],  "plots", "enrichmentAtRecall", "pairwiseComparisons.Recall{recall}.GTExTissue{GTExTissue}.tsv")
	resources:
		mem_mb = determine_mem_mb
	conda:
		os.path.join(config["envDir"], "eQTLEnv.yml")
	script:
		os.path.join(config["codeDir"], "visualization", "plot_enrichment_with_ci_at_recall.R")

# generate final comparison plot (triple boxplot)
enrichmentTables_distance=  [os.path.join(config["outDir"], method, "enrichmentTables", f"enrichmentTable.{distance_min}to{distance_max}Kb.tsv") for distance_min, distance_max in zip(config["distances_min"], config["distances_max"]) for method in config['methods']]
predTables_distance = [os.path.join(config["outDir"], method, "recallTables", "recallTable.byDistance.tsv") for method in config['methods']]
maps = [os.path.join(config["outDir"], method, "intermediate", "GTExTissueBiosampleMap.tsv") for method in config["methods"]]

rule plot_thresholded_performance_comparison:
	input:
		enrichmentTable_files = enrichmentTables_distance,
		predTable_files = predTables_distance,
		map_files = maps,
		colorPalette = os.path.join(config["outDir"], "plots", "colorPalette.tsv"),
	params:
		distances_min = config["distances_min"],
		distances_max = config["distances_max"]
	output:
		outFile = os.path.join(config["outDir"], "plots", "thresholdedPerformanceComparison.pdf"),
		enrAllTable = os.path.join(config["outDir"], "plots", "allMatchedEnrichments.tsv"),
		predictionMetrics = os.path.join(config["outDir"], "plots", "allMatchedPredictionMetrics.tsv"),
		outScatter = os.path.join(config["outDir"], "plots", "allMatchedThresholdedPerformance.scatter.pdf")
	resources:
		mem_mb = determine_mem_mb
	conda:
		os.path.join(config["envDir"], "eQTLEnv.yml")
	script:
		os.path.join(config["codeDir"], "visualization", "plot_thresholded_performance_comparison.R")

# plot all thresholded enrichments + recalls
rule plot_all_thresholded_metrics:
	input:
		enrichmentTable_files = enrichmentTables_distance,
		predTable_files = predTables_distance,
		colorPalette = os.path.join(config["outDir"], "plots", "colorPalette.tsv"),
	params:
		distances_min = config["distances_min"],
		distances_max = config["distances_max"]
	output:
		outFile = os.path.join(config["outDir"], "plots", "allThresholdedEnrichmentRecall.byDistance.pdf"),
	resources:
		mem_mb = determine_mem_mb
	conda:
		os.path.join(config["envDir"], "eQTLEnv.yml")
	script:
		os.path.join(config["codeDir"], "visualization", "plot_all_thresholded_enrichment_recall.R")

# plot n variants per tissue
rule plot_variants_per_tissue:
	input:
		variantsPerTissue =  [os.path.join(config["outDir"], method, "intermediate", "nVariantsPerGTExTissue.tsv") for method in config["methods"]]
	params:
		distances_min = config["distances_min"],
		distances_max = config["distances_max"],
		methods = config["methods"]
	output:
		out_plot = os.path.join(config["outDir"], "plots", "variantsPerTissue.pdf"),
		out_table = os.path.join(config["outDir"], "plots", "avgVariantsPerTissue.tsv")
	resources:
		mem_mb = determine_mem_mb
	conda:
		os.path.join(config["envDir"], "eQTLEnv.yml")
	script:
		os.path.join(config["codeDir"], "visualization", "plot_number_variants.R")

# heatmaps
rule plot_enrichment_heatmaps:
	input:
		enrichmentTable = os.path.join(config["outDir"], "{method}", "enrichmentTables", "enrichmentTable.0to30000Kb.tsv"),
		enhancerSizes =  os.path.join(config["outDir"], "{method}", "intermediate", "basesPerEnhancerSet.tsv")
	params:
		p_threshold = config["thresholdPval"]
	output:
		outFile_combined =  os.path.join(config["outDir"], "plots", "enrichmentHeatmaps", "{method}.enrichmentHeatmap.withMetrics.pdf"),
		outFile_alone = os.path.join(config["outDir"], "plots", "enrichmentHeatmaps", "{method}.enrichmentHeatmap.pdf")
	resources:
		mem_mb = determine_mem_mb
	conda:
		os.path.join(config["envDir"], "eQTLEnv.yml")
	script:
		os.path.join(config["codeDir"], "visualization", "plot_enrichment_heatmap.R")

rule plot_recall_heatmaps:
	input:
		recallTable = os.path.join(config["outDir"], "{method}", "recallTables", "recallTable.byDistance.tsv"),
		enhancerSizes =  os.path.join(config["outDir"], "{method}", "intermediate", "basesPerEnhancerSet.tsv")
	output:
		outFile_combined =  os.path.join(config["outDir"], "plots", "recallHeatmaps", "{method}.recallHeatmap.withMetrics.pdf"),
		outFile_alone = os.path.join(config["outDir"], "plots", "recallHeatmaps", "{method}.recallHeatmap.pdf")
	resources:
		mem_mb = determine_mem_mb
	conda:
		os.path.join(config["envDir"], "eQTLEnv.yml")
	script:
		os.path.join(config["codeDir"], "visualization", "plot_recall_heatmap.R")


# html report
rule generate_html_report:
	input:
		colorPalette = os.path.join(config["outDir"], "plots", "colorPalette.tsv"),
		enrAllTable = os.path.join(config["outDir"], "plots", "allMatchedEnrichments.tsv"),
		predictionMetrics = os.path.join(config["outDir"], "plots", "allMatchedPredictionMetrics.tsv"),
		er_combined = expand(os.path.join(config["outDir"],  "plots", "enrichmentRecall", "enrichmentRecall.GTExTissue{GTExTissue}.tsv"), GTExTissue=GTExTissues_plus_all),
		enrMatrices_CRISPRthresh = expand(os.path.join(config["outDir"], "{method}", "enrichmentTables", "enrichmentTable.0to30000Kb.tsv"), method=config["methods"]),
		nVariants = os.path.join(config["outDir"], "plots", "avgVariantsPerTissue.tsv")
	params:
		methods = config["methods"],
		distances_min = config["distances_min"],
		distances_max = config["distances_max"],
		distances = config["distances"],
		GTExTissues_matched = GTExTissues_plus_all
	resources:
		mem_mb = determine_mem_mb
	output:
		htmlReport = os.path.join(config["outDir"], "benchmarkingReport.html")
	conda:
		os.path.join(config["envDir"], "eQTLEnv.yml")
	script:
		os.path.join(config["codeDir"], "visualization", "benchmarking_report.Rmd")






