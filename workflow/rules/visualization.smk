# color palette 
# generate color palette (run once overall)
colors = []
for x in config["methods"]:
	# list of color inputs
	colors.append(methods_config.loc[x, "color"])
	
rule generate_color_palette:
	params:
		user_inputs = colors,
		names = config["methods"],
		methods_config = config["methodsTable"]
	output:
		colorPalette = os.path.join(config["outDir"], "colorPalette.rds"),
		colorPalettePlotting = os.path.join(config["outDir"], "colorPalettePlotting.rds")
		
	conda: 
		os.path.join(config["envDir"], "eQTLEnv.yml")
	script: 
		os.path.join(config["codeDir"], "color_palette.R")

# gather data for enrichment recall curve
# for each GTEx tissue with at least one matched biosample
rule gather_enrichment_recall:
	input: 
		predTable = os.path.join(config["outDir"], "{method}", "predictionTables", "GTExTissue{GTExTissue}.Biosample{biosample}.byThreshold.tsv"),
		enrichmentTable = os.path.join(config["outDir"], "{method}", "enrichmentTables", "giant_enrichmentTable.threshold.tsv"),
		thresholdSpan = os.path.join(config["outDir"], "{method}", "thresholdSpan.tsv")

	params:
		outDir = config["outDir"],
	conda: 
		os.path.join(config["envDir"], "eQTLEnv.yml")
	output:
		ERCurveTable = os.path.join(config["outDir"], "{method}", "enrichmentRecallTables", "GTExTissue{GTExTissue}.Biosample{biosample}.tsv")
	script:
		os.path.join(config["codeDir"], "gather_enrichment_recall.R")


# combined enrichment recall curves-- one per GTEx tissue-- across methods
rule plot_enrichment_recall_curve:
	input:
		colorPalette = os.path.join(config["outDir"], "colorPalettePlotting.rds"),
		# read in ER tables within script using methods_config
	params:
		codeDir = config["codeDir"],
		outDir = config["outDir"],
		methods = config["methods"],
		methods_config = config["methodsTable"]
	output:
		er_combined = os.path.join(config["outDir"],  "plots", "enrichmentRecall.GTExTissue{GTExTissue}.pdf"),
		er_combined_table = os.path.join(config["outDir"],  "plots", "enrichmentRecall.GTExTissue{GTExTissue}.tsv")
	conda:
		os.path.join(config["envDir"], "eQTLEnv.yml")
	script:
		os.path.join(config["codeDir"], "combined_er_curves_integrated.R")
		
# plot enrichment with confidence intervals across methods, once per GTEx tissue and user-defined recall
rule plot_enrichment_with_ci:
	input:
		colorPalette = os.path.join(config["outDir"], "colorPalettePlotting.rds"),
		# read in ER tables within script using methods_config
	params:
		codeDir = config["codeDir"],
		outDir = config["outDir"],
		methods = config["methods"],
		methods_config = config["methodsTable"],
	output:
		enr_at_recall = os.path.join(config["outDir"],  "plots", "enrichments.Recall{recall}.GTExTissue{GTExTissue}.pdf"),
		enr_at_recall_table = os.path.join(config["outDir"],  "plots", "enrichments.Recall{recall}.GTExTissue{GTExTissue}.tsv"),
		sign_table = os.path.join(config["outDir"],  "plots", "pairwiseComparisons.Recall{recall}.GTExTissue{GTExTissue}.tsv")
	conda:
		os.path.join(config["envDir"], "eQTLEnv.yml")
	script:
		os.path.join(config["codeDir"], "enrichment_with_ci_at_recall.R")


# generate final comparison plot (triple boxplot)
rule plot_final_comparison_grouped_by_distance:
	input:
		enrichmentTables = expand(os.path.join(config["outDir"], "{method}", "enrichmentTables", "enrichmentTable.under{distance}bp.tsv"), method=config["methods"], distance=config["distances"]),
		# read in prediction tables individually within method using biosample/tissue pairings
		colorPalette = os.path.join(config["outDir"], "colorPalettePlotting.rds"),
	params:
		methods = config["methods"],
		distances = config["distances"],
		outDir = config["outDir"],
		methods_config = config["methodsTable"]
		
	output:
		outFile = os.path.join(config["outDir"], "plots", "final_comparison_figure_grouped_by_dist.pdf"),
		enrAllTable = os.path.join(config["outDir"], "plots", "all_matched_enrichments.tsv"),
		predictionMetrics = os.path.join(config["outDir"], "plots", "all_matched_prediction_metrics.tsv")
	conda:
		os.path.join(config["envDir"], "eQTLEnv.yml")
	script:
		os.path.join(config["codeDir"], "combined_plots_grouped_by_dist.R")

# plot averge enrichment/recall across tissue/biosample matches
rule aggregate_enrichment_recall:
	input:
		# all enrichment-recall tables for each GTEx tissue
		enrichmentRecallCurves = expand(os.path.join(config["outDir"],  "plots", "enrichmentRecall.GTExTissue{GTExTissue}.tsv"), GTExTissue=GTExTissues_matched),
		colorPalette = os.path.join(config["outDir"], "colorPalettePlotting.rds"),
	params:
		methods = config["methods"],
		outDir = config["outDir"],
		methods_config = config["methodsTable"],
		matchedGTExTissues = GTExTissues_matched
	output:
		er_aggregate_plot = os.path.join(config["outDir"], "plots", "aggregate_enrichment_recall.pdf"),
		er_aggregate_table = os.path.join(config["outDir"], "plots", "aggregate_enrichment_recall.tsv"),
	conda:
		os.path.join(config["envDir"], "eQTLEnv.yml")
	script:
		os.path.join(config["codeDir"], "aggregate_enrichment_recall.R")

# html report
rule generate_html_report:
	input:
		colorPalette = os.path.join(config["outDir"], "colorPalettePlotting.rds"),
		enrAllTable = os.path.join(config["outDir"], "plots", "all_matched_enrichments.tsv"),
		predictionMetrics = os.path.join(config["outDir"], "plots", "all_matched_prediction_metrics.tsv"),
		er_combined = expand(os.path.join(config["outDir"],  "plots", "enrichmentRecall.GTExTissue{GTExTissue}.tsv"), GTExTissue=GTExTissues_matched),
		enrMatrices_CRISPRthresh = expand(os.path.join(config["outDir"], "{method}", "enrichmentTables", "enrichmentTable.under30000000bp.tsv"), method=config["methods"])
	params:
		methods = config["methods"],
		distances = config["distances"],
		outDir = config["outDir"],
		methods_config = config["methodsTable"],
		matchedGTExTissues = GTExTissues_matched
	output:
		htmlReport = os.path.join(config["outDir"], "benchmarking_report.html")
	conda:
		os.path.join(config["envDir"], "eQTLEnv.yml")
	script:
		os.path.join(config["codeDir"], "benchmarking_report.Rmd")


#############################################################




