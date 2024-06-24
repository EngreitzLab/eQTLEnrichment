# return file for variant counts per tissue for each distance stratification
# variant list columns: 1-3 (loc), 4 (variantID), 5 (gene), 6 (tissue), 7 (PIP),  8 (distance group)
rule get_variants_per_GTEx_tissue:
	input:
		filteredGTExVariantsFinal = os.path.join(config["outDir"], "{method}", "intermediate", "GTExVariants.filteredForMethod.tsv.gz")
	params:
		distances_max= config["distances_max"]
	output: 
		variantsPerTissue =  os.path.join(config["outDir"], "{method}", "intermediate", "nVariantsPerGTExTissue.tsv"),
	resources:
		mem_mb = determine_mem_mb
	conda: 
		os.path.join(config["envDir"], "eQTLEnv.yml")
	script: 
		os.path.join(config["codeDir"], "enrichment_recall", "get_variants_per_tissue.R")

# compute count overlap of variants and predictions with variants stratified by distance and predictions thresholded as input
rule compute_count_matrix_by_distance:
	input:
		variantsPerTissue =  os.path.join(config["outDir"], "{method}", "intermediate", "nVariantsPerGTExTissue.tsv"),
		varPredInt = lambda wildcards: [os.path.join(config["outDir"], wildcards.method, "biosamples", biosample, "GTExVariants-enhancerPredictionsInt.tsv.gz") for biosample in methods_config.loc[wildcards.method, "biosamples"]]
	params:
		threshold = lambda wildcards: methods_config.loc[wildcards.method, "threshold"],
		biosamples = lambda wildcards: methods_config.loc[wildcards.method, "biosamples"]
	output:
		countMatrix = temp(os.path.join(config["outDir"], "{method}", "countMatrices", "count_matrix.{distance_min}to{distance_max}Kb.tsv"))
	resources:
		mem_mb = determine_mem_mb
	conda:
		os.path.join(config["envDir"], "eQTLEnv.yml")
	script:
		os.path.join(config["codeDir"], "enrichment_recall", "count_matrix_distance.R")

# compute count overlap across a range of prediction thresholds
rule compute_count_matrix_by_threshold:
	input:
		thresholdSpan = os.path.join(config["outDir"], "{method}", "intermediate", "thresholdSpan.tsv"),
		variantsPerTissue =  os.path.join(config["outDir"], "{method}", "intermediate", "nVariantsPerGTExTissue.tsv"),
		varPredInt = lambda wildcards: [os.path.join(config["outDir"], wildcards.method, "biosamples", biosample, "GTExVariants-enhancerPredictionsInt.tsv.gz") for biosample in methods_config.loc[wildcards.method, "biosamples"]],
		map = os.path.join(config["outDir"], "{method}", "intermediate", "GTExTissueBiosampleMap.tsv")
	params:
		biosamples = lambda wildcards: methods_config.loc[wildcards.method, "biosamples"],
	output:
		countMatrix = temp(os.path.join(config["outDir"], "{method}", "countMatrices", "count_matrix.acrossThresholds.tsv"))
	resources:
		mem_mb = determine_mem_mb
	conda:
		os.path.join(config["envDir"], "eQTLEnv.yml")
	script:
		os.path.join(config["codeDir"], "enrichment_recall", "count_matrix_across_thresholds.R")

# generate matrix with enrichment values for each GTEx tissue/biosample intersection for each distance
rule compute_enrichment_matrix_by_distance:
	input: 
		countMatrix = os.path.join(config["outDir"], "{method}", "countMatrices", "count_matrix.{distance_min}to{distance_max}Kb.tsv"),
		variantsPerGTExTissueByDist = os.path.join(config["outDir"], "{method}", "intermediate", "nVariantsPerGTExTissue.tsv"),
		commonVarCount = os.path.join(config["outDir"], "variants", "distalNoncodingBackgroundSNPCount.txt"),
		commonVarInt = lambda wildcards: [os.path.join(config["outDir"], wildcards.method, "biosamples", biosample, "distalNoncodingBackgroundSNPs-enhancerPredictionsInt.tsv.gz") for biosample in methods_config.loc[wildcards.method, "biosamples"]]
	params:
		biosamples = lambda wildcards: methods_config.loc[wildcards.method, "biosamples"],
		threshold = lambda wildcards: methods_config.loc[wildcards.method, "threshold"],
		distances_max = config["distances_max"],
		thresholdPval = config["thresholdPval"]
	output: 
		enrichmentTable = os.path.join(config["outDir"], "{method}", "enrichmentTables", "enrichmentTable.{distance_min}to{distance_max}Kb.tsv")
	resources:
		mem_mb = determine_mem_mb
	conda: 
		os.path.join(config["envDir"], "eQTLEnv.yml")
	script: 
		os.path.join(config["codeDir"], "enrichment_recall", "counts_to_enrichment_distance.R")

# generate matrix with enrichment values for each GTEx tissue/biosample intersection for each threshold
rule compute_enrichment_matrix_by_threshold:
	input: 
		countMatrix = os.path.join(config["outDir"], "{method}", "countMatrices", "count_matrix.acrossThresholds.tsv"),
		variantsPerGTExTissue = os.path.join(config["outDir"], "{method}", "intermediate", "nVariantsPerGTExTissue.tsv"),
		thresholdSpan = os.path.join(config["outDir"], "{method}", "intermediate", "thresholdSpan.tsv"),
		commonVarCount = os.path.join(config["outDir"], "variants", "distalNoncodingBackgroundSNPCount.txt"),
		map = os.path.join(config["outDir"], "{method}", "intermediate", "GTExTissueBiosampleMap.tsv"),
		commonVarInt = lambda wildcards: [os.path.join(config["outDir"], wildcards.method, "biosamples", biosample, "distalNoncodingBackgroundSNPs-enhancerPredictionsInt.tsv.gz") for biosample in methods_config.loc[wildcards.method, "biosamples"]]
	params:
		biosamples = lambda wildcards: methods_config.loc[wildcards.method, "biosamples"],
		thresholdPval = config["thresholdPval"]
	output: 
		enrichmentTable = os.path.join(config["outDir"], "{method}", "enrichmentTables", "enrichmentTable.acrossThresholds.tsv.gz")
	resources:
		mem_mb = determine_mem_mb
	conda: 
		os.path.join(config["envDir"], "eQTLEnv.yml")
	script: 
		os.path.join(config["codeDir"], "enrichment_recall", "counts_to_enrichment_across_thresholds.R")


# compute number of base pairs in enhancers above provided
rule compute_enhancer_set_size:
	input:
		predictionsSorted = lambda wildcards: [os.path.join(config["outDir"], wildcards.method, "biosamples",  biosample, "enhancerPredictions.sorted.bed.gz") for biosample in methods_config.loc[wildcards.method, "biosamples"]],
	params:
		biosamples = lambda wildcards: methods_config.loc[wildcards.method, "biosamples"],
		threshold = lambda wildcards: methods_config.loc[wildcards.method, "threshold"]
	output:
		basesPerEnhancerSet = os.path.join(config["outDir"], "{method}", "intermediate", "basesPerEnhancerSet.tsv")
	resources:
		mem_mb = determine_mem_mb
	conda: 
		os.path.join(config["envDir"], "eQTLEnv.yml")
	shell:
		"""			
		set +o pipefail;

		# Create an output file with a header
		echo -e "biosample\tbp" > {output.basesPerEnhancerSet}

		# Convert biosamples and predFiles to arrays
		IFS=' ' read -r -a biosampleArray <<< "{params.biosamples}"
		IFS=' ' read -r -a predFileArray <<< "{input.predictionsSorted}"

		# Loop over the prediction files and corresponding biosamples
		for i in "${{!predFileArray[@]}}"
		do
			biosample=${{biosampleArray[$i]}}
			pred=${{predFileArray[$i]}}

			# Calculate the metric for the current prediction file
			metric=$(zcat $pred | awk -v threshold={params.threshold} '$6 >= threshold' | cut -f 1-3 | bedtools merge -i stdin | awk 'BEGIN {{FS=OFS="\t"}} {{print $3-$2}}' | awk '{{s+=$1}} END {{print s}}')

			# Append the biosample and metric to the output file
			echo -e "${{biosample}}\t${{metric}}" >> {output.basesPerEnhancerSet}
		done
		"""