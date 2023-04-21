# return files for all variants and for each distance stratification
# variant list columns: 1-3 (loc), 4 (hgID), 5 (tissue), 6 (ens_id), 7 (PIP), 8 (TPM), 9 (distance)
rule get_variants_per_GTEx_tissue:
	input:
		filteredGTExVariantsFinal = os.path.join(config["outDir"], "{method}", "GTExVariants.filteredForMethod.tsv")
	params:
		distances = config["distances"],
		GTExTissues=config["GTExTissues"],
		variantsPerGTExTissueByDist = expand(os.path.join(config["outDir"], "{{method}}", "nVariantsPerGTExTissue.under{distance}bp.tsv"), distance=config["distances"])
	output: 
		variantsPerGTExTissue = os.path.join(config["outDir"], "{method}", "nVariantsPerGTExTissue.tsv"),
	conda: 
		os.path.join(config["envDir"], "eQTLEnv.yml")
	script:
		os.path.join(config["codeDir"], "get_variants_per_tissue.R")

rule compute_count_matrix_by_distance:
	input:
		threshold = lambda wildcards: methods_config.loc[wildcards.method, "threshold"]
		# read in wthin script: variantsPredictionsInt = os.path.join(config["outDir"], "{method}", "{biosample}", "GTExVariants-enhancerPredictionsInt.tsv.gz"),
	params:
		biosamples=lambda wildcards: methods_config.loc[wildcards.method, "biosamples"],
		GTExTissues=config["GTExTissues"],
		outDir = config["outDir"],
	output:
		countMatrix = os.path.join(config["outDir"], "{method}", "count_matrix.under{distance}bp.tsv")
	conda:
		os.path.join(config["envDir"], "eQTLEnv.yml")
	script:
		os.path.join(config["codeDir"], "count_matrix_distance.R")

rule compute_count_matrix_by_threshold:
	input:
		biosamples=lambda wildcards: methods_config.loc[wildcards.method, "biosamples"]
		# read in wthin script: variantsPredictionsInt = os.path.join(config["outDir"], "{method}", "{biosample}", "GTExVariants-enhancerPredictionsInt.tsv.gz"),
	params:
		GTExTissues=config["GTExTissues"],
		outDir = config["outDir"]
	output:
		countMatrix = os.path.join(config["outDir"], "{method}", "count_matrix.thresh{threshold}.tsv")
	conda:
		os.path.join(config["envDir"], "eQTLEnv.yml")
	script:
		os.path.join(config["codeDir"], "count_matrix_threshold.R")

# generate matrix with enrichment values for each GTEx tissue/biosample intersection for each distance
rule compute_enrichment_matrix_by_distance:
	input: 
		countMatrix = os.path.join(config["outDir"], "{method}", "count_matrix.under{distance}bp.tsv"),
		variantsPerGTExTissueByDist = os.path.join(config["outDir"], "{method}", "nVariantsPerGTExTissue.under{distance}bp.tsv"),
		# read in commonVarPredictionsInts files per biosample within Rscript
	params:
		codeDir = config["codeDir"],
		outDir = config["outDir"],
		GTExTissues = config["GTExTissues"],
		biosamples = lambda wildcards: methods_config.loc[wildcards.method, "biosamples"],
		nCommonVariants = config["nBGVariants"],
		threshold = lambda wildcards: methods_config.loc[wildcards.method, "threshold"],
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
		# read in commonVarPredictionsInts files per biosample within Rscript
	params:
		codeDir = config["codeDir"],
		outDir = config["outDir"],
		GTExTissues = config["GTExTissues"],
		biosamples = lambda wildcards: methods_config.loc[wildcards.method, "biosamples"],
		nCommonVariants = config["nBGVariants"]
	output: 
		enrichmentTable = os.path.join(config["outDir"], "{method}", "enrichmentTable.thresh{threshold}.tsv")
	conda: 
		os.path.join(config["envDir"], "eQTLEnv.yml")
	script: 
		os.path.join(config["codeDir"], "counts_to_enrichment_threshold.R")