
rule compute_prediction_table_by_distance:
	input:
		variantsPredictionsInt = os.path.join(config["outDir"], "{method}", "{Biosample}", "GTExVariants-enhancerPredictionsInt.tsv.gz"),
		filteredGTExVariantsFinal = os.path.join(config["outDir"], "{method}", "GTExVariants.filteredForMethod.tsv")
	params:
		distances = config["distances"],
		threshold = lambda wildcards: methods_config.loc[wildcards.method, "threshold"],
	output:
		predTable = os.path.join(config["outDir"], "{method}", predictionTables, "GTExTissue{GTExTissue}.Biosample{Biosample}.byDistance.tsv")
	conda:
			os.path.join(config["envDir"], "eQTLEnv.yml")
	script:
		os.path.join(config["codeDir"], "compute_prediction_table_distance.R")


rule compute_prediction_table_by_threshold:
	input:
		variantsPredictionsInt = os.path.join(config["outDir"], "{method}", "{Biosample}", "GTExVariants-enhancerPredictionsInt.tsv.gz"),
		filteredGTExVariantsFinal = os.path.join(config["outDir"], "{method}", "GTExVariants.filteredForMethod.tsv")
	params:
		thresholds = lambda wildcards: methods_config.loc[wildcards.method, "thresholdSpan"]
	output:
		predTable = os.path.join(config["outDir"], "{method}", predictionTables, "GTExTissue{GTExTissue}.Biosample{Biosample}.byThreshold.tsv")
	conda:
			os.path.join(config["envDir"], "eQTLEnv.yml")
	script:
		os.path.join(config["codeDir"], "compute_prediction_table_threshold.R")



######################### OLD ############################

# load prediction method config file
#methods_config_file = config["methodsTable"]
#methods_config = pd.read_table(methods_config_file, na_values="").fillna("None").set_index("method", drop=False)
#methods_config["GTExTissue_map"] = methods_config["GTExTissue_map"].apply(eval)
#methods_config["biosample_map"] = methods_config["biosample_map"].apply(eval)

# generate output files
variantsByTissueFiles = []
predictionsByBiosampleFiles = []
predTablesFiles = []
thresholdTableFiles = []

#for x in config["methods"]:
#	variantsByTissueFiles.extend(expand(os.path.join(config["outDir"], x, "eGenePrediction",  #"{GTExTissue}.{biosample}.filteredVariants.tsv"), zip, GTExTissue=methods_config.loc[x, "GTExTissue_map"], #biosample=methods_config.loc[x, "biosample_map"]))
#	predictionsByBiosampleFiles.extend(expand(os.path.join(config["outDir"], x, "{biosample}", #"enhancerPredictions.sorted.bed.gz"), biosample=methods_config.loc[x, "biosample_map"]))
#	predTablesFiles.extend(expand(os.path.join(config["outDir"], x, "eGenePrediction", #"{GTExTissue}.{Biosample}.predictionTable.tsv"), zip, GTExTissue=methods_config.loc[x, "GTExTissue_map"], #Biosample=methods_config.loc[x, "biosample_map"]))
#	thresholdTableFiles.extend(expand(os.path.join(config["outDir"], "thresholdTables", x, "{GTExTissue}.{biosample}.tsv"), zip, #GTExTissue=methods_config.loc[x, "GTExTissue_map"], biosample=methods_config.loc[x, "biosample_map"]))

	
# split variants to tissue (just the ones needed for predictions)
rule split_variants_by_tissue:
	input: 
		filteredGTExVariantsFinal = os.path.join(config["outDir"], "{method}", "GTExVariants.filteredForMethod.tsv")
	params:
		codeDir = config["codeDir"],
		outDir = config["outDir"]
	output:
		variantsByTissue = os.path.join(config["outDir"], "{method}", "eGenePrediction", "GTExTissue{GTExTissue}.Biosample{biosample}.filteredVariants.tsv")
	conda: 
		os.path.join(config["envDir"], "eQTLEnv.yml")
	shell:
		"""
		set +o pipefail;
		# input variant file columns: 1-3 (loc), 4 (hgID), 5 (tissue), 6 (gene hgnc), 7 (PIP)
		cat {input.filteredGTExVariantsFinal} | awk -v var={wildcards.GTExTissue} '$5==var' | sort -k1,1 -k2,2n > {output.variantsByTissue}
		"""

# generate prediction table and threshold tables
rule make_prediction_table_for_calcs:
	input:
		variantsByTissue = os.path.join(config["outDir"], "{method}", "eGenePrediction", "GTExTissue{GTExTissue}.Biosample{biosample}.filteredVariants.tsv"),
		predictionsByBiosample =  os.path.join(config["outDir"], "{method}", "{biosample}", "enhancerPredictions.sorted.bed.gz"),
		
	params:
		outDir = config["outDir"],
		codeDir = config["codeDir"],
		span = lambda wildcards: methods_config.loc[wildcards.method, "thresholdSpan"]
	output:
		predTable = os.path.join(config["outDir"], "{method}", "eGenePrediction", "GTExTissue{GTExTissue}.Biosample{biosample}.predictionTable.forCalcs.tsv"),
		thresholdTable = os.path.join(config["outDir"], "thresholdTables", "{method}", "GTExTissue{GTExTissue}.Biosample{biosample}.tsv")
	conda: 
		os.path.join(config["envDir"], "eQTLEnv.yml")
	shell:
		"""
		set +o pipefail;
		
		# intersect predictions with relevant variants 
		cat {input.variantsByTissue} | cut -f 1-4 | sed 1d | sort -k 1,1 -k2,2n | uniq > {params.outDir}/{wildcards.method}/eGenePrediction/{wildcards.GTExTissue}.temp.uniqueVariants.tsv		
		
		zcat {input.predictionsByBiosample} | bedtools intersect -wa -wb -a {params.outDir}/{wildcards.method}/eGenePrediction/{wildcards.GTExTissue}.temp.uniqueVariants.tsv -b stdin | cut -f 4,9,10 > {params.outDir}/{wildcards.method}/eGenePrediction/{wildcards.GTExTissue}.temp.predictionTargetGenes.tsv
		
		# classify variants, generate tables
		Rscript {params.codeDir}/classify_enhancer_predictions_calcs.R --variants {input.variantsByTissue} --pred {params.outDir}/{wildcards.method}/eGenePrediction/{wildcards.GTExTissue}.temp.predictionTargetGenes.tsv --outDir {params.outDir} --outThresh {output.thresholdTable}  --biosample {wildcards.biosample} --span "{params.span}" > {output.predTable}
		
		rm {params.outDir}/{wildcards.method}/eGenePrediction/{wildcards.GTExTissue}.temp.uniqueVariants.tsv
		rm {params.outDir}/{wildcards.method}/eGenePrediction/{wildcards.GTExTissue}.temp.predictionTargetGenes.tsv
		"""


