# threshold predictors by score
rule threshold_predictions:
	input: 
		predictionsSorted = os.path.join(config["outDir"], "{method}", "{biosample}", "enhancerPredictions.sorted.bed.gz"),
	params:
		ourDir = config["outDir"],
	conda: 
		os.path.join(config["envDir"], "eQTLEnv.yml")
	output:
		predictionsThresholded = os.path.join(config["outDir"], "{method}", "{biosample}", "enhancerPredictions.thresholded.{threshold}.bed")
	script:
		os.path.join(config["codeDir"], "threshold_predictions.R")

# gzip thresholds
rule gzip_thresholds:
	input:
		predictionsThresholded = os.path.join(config["outDir"], "{method}", "{biosample}", "enhancerPredictions.thresholded.{threshold}.bed")
	output:
		predictionsThresholdedGzipped = os.path.join(config["outDir"], "{method}", "{biosample}", "enhancerPredictions.thresholded.{threshold}.bed.gz")
	shell:
		"""
			set +o pipefail;
	
			gzip -f {input.predictionsThresholded} > {output.predictionsThresholdedGzipped}
		"""
		
# gather data for and plot enrichment recall curves
rule plot_enrichment_recall:
	input: 
		# threshold table files for a given method and all biosample-tissue pairings
		thresholdTables = lambda wildcards: expand(os.path.join("thresholdTables", wildcards.method, "{GTExTissue}.{biosample}.tsv"), zip, GTExTissue=methods_config.loc["{method}", "GTExTissue_map"], biosample=methods_config.loc["{method}", "biosamples"]),
		# enrichment tables for given method and all threhsolds
		enrichmentTables = lambda wildcards: expand(os.path.join(config["outDir"], wildcards.method, "enrichmentTable.{threshold}.tsv", threshold=methods_config.loc["{method}", "thresholdSpan"]))
	params:
		ourDir = config["outDir"]
	conda: 
		os.path.join(config["envDir"], "eQTLEnv.yml")
	output:
		ERCurveMeanFull = os.path.join(config["outDir"], "{method}", "ERCurveMean.full.pdf"),
		ERCurveMeanZoom = os.path.join(config["outDir"], "{method}", "ERCurveMean.zoom.pdf"),
		ERCurveMaxFull = os.path.join(config["outDir"], "{method}", "ERCurveMax.full.pdf"),
		ERCurveMaxZoom = os.path.join(config["outDir"], "{method}", "ERCurveMax.zoom.pdf"),
		ERCurveTable = os.path.join(config["outDir"], "{method}", "ERCurveTable.tsv")
	script:
		os.path.join(config["codeDir"], "plot_enrichment_recall.R")