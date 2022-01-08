
# threshold baseline predictors by score; otherwise, return unchanged file
rule threshold_predictions:
	input: 
		predictionsSorted = os.path.join(config["outDir"], "{method}", "{biosample}", "enhancerPredictions.sorted.bed.gz"),
		thresholds = os.path.join(config["outDir"], "thresholdTables", "calculatedThresholds.tsv")
	params:
		ourDir = config["outDir"],
		userThresh = lambda wildcards: methods_config.loc[wildcards.method, "threshold"]
	conda: 
		os.path.join(config["envDir"], "eQTLEnv.yml")
	output:
		predictionsThresholded = os.path.join(config["outDir"], "{method}", "{biosample}", "enhancerPredictions.thresholded.bed")
	script:
		os.path.join(config["codeDir"], "threshold_predictions.R")

# gzip thresholds
rule gzip_thresholds:
	input:
		predictionsThresholded = os.path.join(config["outDir"], "{method}", "{biosample}", "enhancerPredictions.thresholded.bed")
	output:
		predictionsThresholdedGzipped = os.path.join(config["outDir"], "{method}", "{biosample}", "enhancerPredictions.thresholded.bed.gz")
	shell:
		"""
			set +o pipefail;
	
			gzip -f {input.predictionsThresholded} > {output.predictionsThresholdedGzipped}
		"""
