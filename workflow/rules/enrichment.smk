# return files for all variants and for each distance stratification
# variant list columns: 1-3 (loc), 4 (hgID), 5 (tissue), 6 (ens_id), 7 (PIP), 8 (TPM), 9 (distance)
rule get_variants_per_GTEx_tissue:
	input:
		filteredGTExVariantsFinal = os.path.join(config["outDir"], "{method}", "GTExVariants.filteredForMethod.tsv")
	params:
		distances = config["distances"],
		GTExTissues=config["GTExTissues"]
	output: 
		variantsPerGTExTissue = os.path.join(config["outDir"], "{method}", "nVariantsPerGTExTissue.tsv"),
		variantsPerGTExTissueByDist = expand(os.path.join(config["outDir"], "{method}", "nVariantsPerGTExTissue.under{distance}bp.tsv"), distances=config["distances"])
	conda: 
		os.path.join(config["envDir"], "eQTLEnv.yml")
	script:
		os.path.join(config["codeDir"], "get_variants_per_tissue.R")

rule compute_count_matrix_by_distance:
	input:
		biosamples=config["biosamples"]
		# read in wthin script: variantsPredictionsInt = os.path.join(config["outDir"], "{method}", "{biosample}", "GTExVariants-enhancerPredictionsInt.tsv.gz"),
	params:
		GTExTissues=config["GTExTissues"],
		outDir = params["outDir"]
	output:
		countMatrix = os.path.join(config["outDir"], "{method}", "count_matrix.under{distance}bp.tsv")
	conda:
		os.path.join(config["envDir"], "eQTLEnv.yml")
	script:
		os.path.join(config["codeDir"], "count_matrix_by_distance.R")

rule compute_count_matrix_by_threshold:
	input:
		biosamples=config["biosamples"]
		# read in wthin script: variantsPredictionsInt = os.path.join(config["outDir"], "{method}", "{biosample}", "GTExVariants-enhancerPredictionsInt.tsv.gz"),
	params:
		GTExTissues=config["GTExTissues"],
		outDir = params["outDir"]
	output:
		countMatrix = os.path.join(config["outDir"], "{method}", "count_matrix.thresh{threshold}.tsv")
	conda:
		os.path.join(config["envDir"], "eQTLEnv.yml")
	script:
		os.path.join(config["codeDir"], "count_matrix_by_distance.R")

# generate matrix with enrichment values for each GTEx tissue/biosample intersection for each distance
rule compute_enrichment_matrix_by_distance:
	input: 
		countMatrix = os.path.join(config["outDir"], "{method}", "count_matrix.under{distance}bp.tsv"),
		variantsPerGTExTissueByDist = os.path.join(config["outDir"], "{method}", "nVariantsPerGTExTissue.under{distance}bp.tsv")
		commonVarPredictionsInt = os.path.join(config["outDir"], "{method}", "{biosample}", "distalNoncodingBackgroundSNPs-enhancerPredictionsInt.tsv.gz")
	params:
		codeDir = config["codeDir"],
		outDir = config["outDir"],
		GTExTissues = config["GTExTissues"],
		biosamples = config["biosamples"],
		nCommonVariants = config["nCommonVariants"]
	output: 
		enrichmentTable = os.path.join(config["outDir"], "{method}", "enrichmentTable.under{distance}bp.tsv")
	conda: 
		os.path.join(config["envDir"], "eQTLEnv.yml")
	script: 
		os.path.join(config["codeDir"], "counts_to_enrichment_distance.R")

# generate matrix with enrichment values for each GTEx tissue/biosample intersection for each threshold
rule compute_enrichment_matrix_by_threshold:
	input: 
		countMatrix = os.path.join(config["outDir"], "{method}", "count_matrix.thresh{threshold}.tsv"),
		variantsPerGTExTissue = os.path.join(config["outDir"], "{method}", "nVariantsPerGTExTissue.tsv"),
		samples = os.path.join(config["outDir"],"{method}", "biosampleList.tsv"),
		commonVarPredictionsInt = os.path.join(config["outDir"], "{method}", "{biosample}", "distalNoncodingBackgroundSNPs-enhancerPredictionsInt.tsv.gz")
	params:
		codeDir = config["codeDir"],
		outDir = config["outDir"],
		GTExTissues = config["GTExTissues"],
		nCommonVariants = config["nCommonVariants"]
	output: 
		enrichmentTable = os.path.join(config["outDir"], "{method}", "enrichmentTable.thresh{threshold}.tsv")
	conda: 
		os.path.join(config["envDir"], "eQTLEnv.yml")
	script: 
		os.path.join(config["codeDir"], "counts_to_enrichment_threshold.R")