
# threshold baseline predictors by score; otherwise, return unchanged file
rule threshold_predictions:
	input: 
		predictionsSorted = os.path.join(config["outDir"], "{method}", "{biosample}", "enhancerPredictions.sorted.bed.gz"),
	params:
		ourDir = config["outDir"],
		threshold = {threshold}
	conda: 
		os.path.join(config["envDir"], "eQTLEnv.yml")
	output:
		predictionsThresholded = os.path.join(config["outDir"], "{method}", "{biosample}", "enhancerPredictions.{threshold}.bed")
	script:
		os.path.join(config["codeDir"], "threshold_predictions.R")

# gzip thresholds
rule gzip_thresholds:
	input:
		predictionsThresholded = os.path.join(config["outDir"], "{method}", "{biosample}", "enhancerPredictions.{threshold}.bed")
	output:
		predictionsThresholdedGzipped = os.path.join(config["outDir"], "{method}", "{biosample}", "enhancerPredictions.{threshold}.bed.gz")
	shell:
		"""
			set +o pipefail;
	
			gzip -f {input.predictionsThresholded} > {output.predictionsThresholdedGzipped}
		"""
