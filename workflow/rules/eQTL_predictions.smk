# load prediction method config file
methods_config_file = config["methodsTable"]
methods_config = pd.read_table(methods_config_file, na_values="").fillna("None").set_index("method", drop=False)
methods_config["GTExTissue_map"] = methods_config["GTExTissue_map"].apply(eval)
methods_config["biosample_map"] = methods_config["biosample_map"].apply(eval)

# generate output files
variantsByTissueFiles = []
predictionsByBiosampleFiles = []
predTablesFiles = []
thresholdTableFiles = []
predThresholdedFiles = []

for x in config["methods"]:
	# list of variant files
	variantsByTissueFiles.extend(expand(os.path.join(config["outDir"], x, "eGenePrediction",  "{GTExTissue}.{biosample}.filteredVariants.tsv"), zip, GTExTissue=methods_config.loc[x, "GTExTissue_map"], biosample=methods_config.loc[x, "biosample_map"]))
	# list of filtered prediction files
	predictionsByBiosampleFiles.extend(expand(os.path.join(config["outDir"], x, "{biosample}", "enhancerPredictions.sorted.bed.gz"), biosample=methods_config.loc[x, "biosample_map"]))
	# list of prediction tables
	predTablesFiles.extend(expand(os.path.join(config["outDir"], x, "eGenePrediction", "{GTExTissue}.{Biosample}.predictionTable.tsv"), zip, GTExTissue=methods_config.loc[x, "GTExTissue_map"], Biosample=methods_config.loc[x, "biosample_map"]))
	# list of variant-prediction intersections
	#predThresholdedFiles.extend(expand(os.path.join(config["outDir"], x, "{biosample}", "enhancerPredictions.{threshold}.bed.gz"), biosample=methods_config.loc[x, "biosamples"], threshold=methods_config.loc[x, "thresholdSpan"]))
	thresholdTableFiles.extend(expand(os.path.join(config["outDir"], "thresholdTables", x, "{GTExTissue}.{biosample}.tsv"), zip, GTExTissue=methods_config.loc[x, "GTExTissue_map"], biosample=methods_config.loc[x, "biosample_map"]))

	
# split variants to tissue (just the ones needed for predictions)
# run once per method x tissue
rule split_variants_by_tissue:
	input: 
		filteredGTExVariantsFinal = os.path.join(config["outDir"], "{method}", "GTExVariants.filteredForMethod.tsv")
	params:
		codeDir = config["codeDir"],
		outDir = config["outDir"]
	output:
		variantsByTissue = os.path.join(config["outDir"], "{method}", "eGenePrediction", "{GTExTissue}.{biosample}.filteredVariants.tsv")
	conda: 
		os.path.join(config["envDir"], "eQTLEnv.yml")
	shell:
		"""
		set +o pipefail;
		# input variant file columns: 1-3 (loc), 4 (hgID), 5 (tissue), 6 (gene hgnc), 7 (PIP)
		cat {input.filteredGTExVariantsFinal} | awk -v var={wildcards.GTExTissue} '$5==var' | sort -k1,1 -k2,2n > {output.variantsByTissue}
		"""
		
# return TSS and gene files to use for proximal mesaures, cell-type specific if provided
rule get_TSS_gene_files:
	input:
		variantsByTissue = os.path.join(config["outDir"], "{method}", "eGenePrediction", "{GTExTissue}.{biosample}.filteredVariants.tsv"),
		geneUniverse = os.path.join(config["outDir"], "{method}", "geneUniverse.tsv")
	params:
		sampleKey = lambda wildcards: methods_config.loc[wildcards.method, "sampleKey"],
		TSS = config["TSS"],
		GTExGeneUniverse = config["GTExGeneUniverse"],
		codeDir = config["codeDir"]
	conda: 
		os.path.join(config["envDir"], "eQTLEnv.yml")
	output:
		specificTSSFile = os.path.join(config["outDir"], "{method}", "eGenePrediction", "{GTExTissue}.{biosample}.TSS.bed"),
		specificGeneFile = os.path.join(config["outDir"], "{method}", "eGenePrediction", "{GTExTissue}.{biosample}.genes.bed")
	script:
		os.path.join(config["codeDir"], "get_TSS_gene_lists.R")
	

# add columns to variant files with proximal genes, proximal TSS, TSS within 100 kb (filter reference gene files to gene universe)
# run once per method x tissue
rule add_proximal_genes:
	input:
		specificTSSFile = os.path.join(config["outDir"], "{method}", "eGenePrediction", "{GTExTissue}.{biosample}.TSS.bed"),
		specificGeneFile = os.path.join(config["outDir"], "{method}", "eGenePrediction", "{GTExTissue}.{biosample}.genes.bed"),
		variantsByTissue = os.path.join(config["outDir"], "{method}", "eGenePrediction", "{GTExTissue}.{biosample}.filteredVariants.tsv"),
	params:
		chrSizes = config["chrSizes"],
		codeDir = config["codeDir"],
		outDir = config["outDir"],	
	output:
		#variantsByTissueProximal = os.path.join(config["outDir"], "{method}", "eGenePrediction", "{GTExTissue}.{biosample}.filteredVariants.proximal.tsv"),
		specificTSSFileSorted = os.path.join(config["outDir"], "{method}", "eGenePrediction", "{GTExTissue}.{biosample}.TSS.bed.sorted"),
		specificGeneFileSorted = os.path.join(config["outDir"], "{method}", "eGenePrediction", "{GTExTissue}.{biosample}.genes.bed.sorted"),
		
	conda: 
		os.path.join(config["envDir"], "eQTLEnv.yml")
	shell:
		"""
		set +o pipefail;
		sort -k 1,1 -k2,2n {input.specificTSSFile} > {output.specificTSSFileSorted}
		sort -k 1,1 -k2,2n {input.specificGeneFile} > {output.specificGeneFileSorted}
		
		# get list of distinct variants, expand variant +/- 100 kb, intersect with location of all TSS; end up with following columns: (1) variant unique ID, (2) gene with TSS within 100 kb
		#cat {input.variantsByTissue} | cut -f 1-4 | sort -k 1,1 -k2,2n | uniq | bedtools slop -i stdin -b 100000 -g {params.chrSizes} | bedtools intersect -wa -wb -a stdin -b {output.specificTSSFileSorted} | cut -f 4,8 > {params.outDir}/{wildcards.method}/eGenePrediction/{wildcards.GTExTissue}.temp.100kbInterval.tsv

		# closest gene body (9), TSS (10), whether eGene has TSS within 100 kb (11)
		# how to handle cases where two genes or TSSs are equidistant? currently reporting only first (can do all, but leads to variants being counted multiple times)
		#bedtools closest -t first -a {input.variantsByTissue} -b {output.specificGeneFileSorted} | cut -f 1-4,6,7,12 | bedtools closest -t first -a stdin -b {output.specificTSSFileSorted} | cut -f 1-7,11 > {params.outDir}/{wildcards.method}/eGenePrediction/{wildcards.GTExTissue}.temp.closest.tsv
		
		#Rscript {params.codeDir}/TSS_close.R --variants {params.outDir}/{wildcards.method}/eGenePrediction/{wildcards.GTExTissue}.temp.closest.tsv --TSSint {params.outDir}/{wildcards.method}/eGenePrediction/{wildcards.GTExTissue}.temp.100kbInterval.tsv > {output.variantsByTissueProximal}
		
		#rm {params.outDir}/{wildcards.method}/eGenePrediction/{wildcards.GTExTissue}.temp.100kbInterval.tsv
		#rm {params.outDir}/{wildcards.method}/eGenePrediction/{wildcards.GTExTissue}.temp.closest.tsv
		
		"""

rule make_prediction_table_for_calcs:
	input:
		variantsByTissue = os.path.join(config["outDir"], "{method}", "eGenePrediction", "{GTExTissue}.{biosample}.filteredVariants.tsv"),
		predictionsByBiosample =  os.path.join(config["outDir"], "{method}", "{biosample}", "enhancerPredictions.sorted.bed.gz"),
		
	params:
		outDir = config["outDir"],
		codeDir = config["codeDir"],
	output:
		predTable = os.path.join(config["outDir"], "{method}", "eGenePrediction", "{GTExTissue}.{biosample}.predictionTable.tsv"),
		thresholdTable = os.path.join(config["outDir"], "thresholdTables", "{method}", "{GTExTissue}.{biosample}.tsv")
	conda: 
		os.path.join(config["envDir"], "eQTLEnv.yml")
	shell:
		"""
		set +o pipefail;
		
		# intersect predictions with relevant variants 
		cat {input.variantsByTissue} | cut -f 1-4 | sed 1d | sort -k 1,1 -k2,2n | uniq > {params.outDir}/{wildcards.method}/eGenePrediction/{wildcards.GTExTissue}.temp.uniqueVariants.tsv		
		
		zcat {input.predictionsByBiosample} | bedtools intersect -wa -wb -a {params.outDir}/{wildcards.method}/eGenePrediction/{wildcards.GTExTissue}.temp.uniqueVariants.tsv -b stdin | cut -f 4,9,10 > {params.outDir}/{wildcards.method}/eGenePrediction/{wildcards.GTExTissue}.temp.predictionTargetGenes.tsv
		
		# classify variants	
		
		Rscript {params.codeDir}/classify_enhancer_predictions_calcs.R --variants {input.variantsByTissue} --pred {params.outDir}/{wildcards.method}/eGenePrediction/{wildcards.GTExTissue}.temp.predictionTargetGenes.tsv --outDir {params.outDir} --outThresh {output.thresholdTable}  --biosample {wildcards.biosample} > {output.predTable}
		
		rm {params.outDir}/{wildcards.method}/eGenePrediction/{wildcards.GTExTissue}.temp.uniqueVariants.tsv
		rm {params.outDir}/{wildcards.method}/eGenePrediction/{wildcards.GTExTissue}.temp.predictionTargetGenes.tsv
		"""


