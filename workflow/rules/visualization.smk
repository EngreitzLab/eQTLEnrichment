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
enrMatrices_gather = []
for x in config["methods"]:
	enrMatrices_gather.extend(expand(os.path.join(config["outDir"], x, "enrichmentTables", "enrichmentTable.thresh{threshold}.tsv"), threshold=methods_config.loc[x, "thresholdSpan"]))
rule gather_enrichment_recall:
	input: 
		predTable = os.path.join(config["outDir"], "{method}", "predictionTables", "GTExTissue{GTExTissue}.Biosample{biosample}.byThreshold.tsv"),
		enrichmentMatrices = enrMatrices_gather
	params:
		outDir = config["outDir"],
		thresholdSpan = lambda wildcards: methods_config.loc[wildcards.method, "thresholdSpan"]
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


# html report
# find set of all GTEx tissues with biosample matches across methods 
GTExTissues_all = []
for x in config["methods"]:
	GTExTissues_all.extend(methods_config.loc[x, "GTExTissue_map"])
GTExTissues_matched = [*set(GTExTissues_all)]

rule generate_html_report:
	input:
		colorPalette = os.path.join(config["outDir"], "colorPalettePlotting.rds"),
		enrAllTable = os.path.join(config["outDir"], "plots", "all_matched_enrichments.tsv"),
		predictionMetrics = os.path.join(config["outDir"], "plots", "all_matched_prediction_metrics.tsv"),
		er_combined = expand(os.path.join(config["outDir"],  "plots", "enrichmentRecall.GTExTissue{GTExTissue}.tsv"), GTExTissue=GTExTissues_matched),
		enrMatrices_CRISPRthresh = expand(os.path.join(config["outDir"], "{method}", "enrichmentTable.under30000000bp.tsv"), method=config["methods"])
	params:
		methods = config["methods"],
		distances = config["distances"],
		outDir = config["outDir"],
		methods_config = methods_config
	output:
		htmlReport = os.path.join(config["outDir"], "benchmarking_report.html")
	conda:
		os.path.join(config["envDir"], "eQTLEnv.yml")
	script:
		os.path.join(config["codeDir"], "benchmarking_report.Rmd")


#############################################################




