# get list of genes considered by GTEx and by the prediction method
rule make_gene_universes:
	input: 
 		methodGeneUniverse = lambda wildcards: methods_config.loc[wildcards.method, "geneUniverse"],
		GTExGeneUniverse = config["GTExGeneUniverse"]
	params:
		outDir = config["outDir"],
		codeDir = config["codeDir"]
	output:
		geneUniverse = os.path.join(config["outDir"], "{method}", "geneUniverse.tsv")
	conda: 
		os.path.join(config["envDir"], "eQTLEnv.yml")
	shell:
			"""
			set +o pipefail;
			
			# filter GTEx genes to method genes
			Rscript {params.codeDir}/process_file_in_various_ways.R --input {input.GTExGeneUniverse} --id_col 4 --genes {input.methodGeneUniverse} --id hgnc --invert False --score_col 6 > {output.geneUniverse}
			"""

# sort enhancer predictions by chromosome & start location & filter to gene universe, get list of samples
rule process_predictions:
	input:
		predFile = lambda wildcards: methods_config.loc[wildcards.method, "predFiles"][wildcards.biosample],
		geneUniverse = os.path.join(config["outDir"], "{method}", "geneUniverse.tsv")
	params:
		codeDir = config["codeDir"],
		outDir = config["outDir"],
		chrSizes = config["chrSizes"],
		scoreCol = lambda wildcards: methods_config.loc[wildcards.method, "score_col"],
		inversePred = lambda wildcards: methods_config.loc[wildcards.method, "inverse_predictor"]
	output:
		predictionsSorted = os.path.join(config["outDir"], "{method}", "{biosample}", "enhancerPredictions.sorted.bed.gz"),
	conda: 
		os.path.join(config["envDir"], "eQTLEnv.yml")
	shell:
		"""
		set +o pipefail;
			
		# sort predictions file
        if [[ {input.predFile} == *.gz ]]
        then
		    zcat {input.predFile} | csvtk cut -t -f chr,start,end,CellType,TargetGene,{params.scoreCol} | sed 1d | bedtools sort -i stdin -faidx {params.chrSizes} > {params.outDir}/{wildcards.method}/{wildcards.biosample}/temp.sortedPred.tsv
		else
            cat {input.predFile} | csvtk cut -t -f chr,start,end,CellType,TargetGene,{params.scoreCol} | sed 1d | bedtools sort -i stdin -faidx {params.chrSizes} > {params.outDir}/{wildcards.method}/{wildcards.biosample}/temp.sortedPred.tsv
        fi

		# invert score if inverted predictor and filter to gene universe and set biosampe column
		Rscript {params.codeDir}/process_file_in_various_ways.R --input {params.outDir}/{wildcards.method}/{wildcards.biosample}/temp.sortedPred.tsv --id_col 5 --genes {input.geneUniverse} --id hgnc --biosample {wildcards.biosample} --biosample_col 4 --invert {params.inversePred} --score_col 6 | gzip > {output.predictionsSorted}

		rm {params.outDir}/{wildcards.method}/{wildcards.biosample}/temp.sortedPred.tsv
			
		"""

# filter GTEx variants by PIP, credible set, and to distal noncoding genes; convert ensembl id to hgnc
rule filter_all_variants:
	input:
		GTExVariants = config["GTExVariants"],
		commonVar = config["bgVariants"],
		partition = config["partition"]
	params: 
		codeDir = config["codeDir"],
		outDir = config["outDir"],
		chrSizes = config["chrSizes"],
		thresholdPIP = config["thresholdPIP"]
	output: 
		filteredGTExVar = os.path.join(config["outDir"], "generalReference", "GTExVariants.filtered.PIPfilt.distalNoncoding.tsv.gz"),
		partitionDistalNoncoding = os.path.join(config["outDir"], "generalReference", "Partition.distalNoncoding.bed"),
		commonVarDistalNoncoding = os.path.join(config["outDir"], "generalReference", "distalNoncodingBackgroundSNPs.bed.gz")
	conda: 
		os.path.join(config["envDir"], "eQTLEnv.yml")	
	shell:
			"""
			set +o pipefail;
			
			# filter partition to distal noncoding
			awk '$4=="ABC" || $4=="AllPeaks" || $4=="Other" || $4=="OtherIntron"' {input.partition} | bedtools sort -i stdin -faidx {params.chrSizes} > {output.partitionDistalNoncoding}
			
			# filter GTEx variants by SUSIE, credible set, PIP; print set of columns: 1-3 (loc), 4 (hgID), 5 (tissue), 6 (ens_id), 7 (PIP); then filter to distal noncoding
			zcat {input.GTExVariants} | awk '$16>={params.thresholdPIP} && $17 != -1  && $9 == "SUSIE"' | cut -f1-4,10,11,16 | bedtools sort -i stdin -faidx {params.chrSizes} | uniq | bedtools intersect -wa -sorted -a stdin -b {output.partitionDistalNoncoding} -g {params.chrSizes} | gzip > {output.filteredGTExVar}

			# filter common variants to distal noncoding
			cat {input.commonVar} | bedtools sort -i stdin -faidx {params.chrSizes} | bedtools intersect -wa -sorted -a stdin -b {output.partitionDistalNoncoding} -g {params.chrSizes}| gzip > {output.commonVarDistalNoncoding}
			"""

# filter variants by expression (based on what tissue they are in and GTEx expression data)
rule filter_variants_by_expression:
	input: 
		exprData = config["GTExExpression"],
		filteredGTExVar = os.path.join(config["outDir"], "generalReference", "GTExVariants.filtered.PIPfilt.distalNoncoding.tsv.gz"),
	params:
		outDir = config["outDir"],
		thresholdTPM = config["thresholdTPM"]
	output:
		GTExExpressedGenes = os.path.join(config["outDir"], "generalReference", "GTExVariants.filtered.PIPfilt.distalNoncoding.expressed.tsv")
	conda: 
		os.path.join(config["envDir"], "eQTLEnv.yml")
	script: 
		os.path.join(config["codeDir"], "GTEx_expression.R")

# add eVariant - eGene TSS variant distance to variant file
# columns: 1-3 (loc), 4 (hgID), 5 (tissue), 6 (ens_id), 7 (PIP), 8 (TPM), 9 (distance group)
rule add_distance_to_variants:
	input:
		GTExExpressedGenes = os.path.join(config["outDir"], "generalReference", "GTExVariants.filtered.PIPfilt.distalNoncoding.expressed.tsv")
	params:
		codeDir = config["codeDir"],
		outDir = config["outDir"],
		TSS = config['TSS'],
		distances = config["distances"]
	output:
		GTExExpressedGenesWithDistance = os.path.join(config["outDir"], "generalReference", "GTExVariants.filtered.PIPfilt.distalNoncoding.expressed.withDistance.tsv")
	conda:
		os.path.join(config["envDir"], "eQTLEnv.yml")
	script:
		os.path.join(config["codeDir"], "add_distance_to_variants.R")


# filter variants to by gene universe for each method
# columns: 1-3 (loc), 4 (hgID), 5 (tissue), 6 (ens_id), 7 (PIP), 8 (TPM), 9 (distance)
rule filter_variants_to_gene_universe:
	input:
		GTExExpressedGenesWithDistance = os.path.join(config["outDir"], "generalReference", "GTExVariants.filtered.PIPfilt.distalNoncoding.expressed.withDistance.tsv"),
		geneUniverse = os.path.join(config["outDir"], "{method}", "geneUniverse.tsv")
	params:
		codeDir = config["codeDir"],
		outDir = config["outDir"]
	output:
		filteredGTExVariantsFinal = os.path.join(config["outDir"], "{method}", "GTExVariants.filteredForMethod.tsv")
	conda: 
		os.path.join(config["envDir"], "eQTLEnv.yml")
	shell:
		"""
			## filter this method's gene universe
			Rscript {params.codeDir}/process_file_in_various_ways.R --input {input.GTExExpressedGenesWithDistance} --id_col 6 --genes {input.geneUniverse} --id hgnc --invert False --score_col 0 > {output.filteredGTExVariantsFinal}
			
		"""

# intersect predictions for each threshold with GTEx variants
# output columns: 1-3 (loc), 4 (hgID), 5 (tissue), 6 (ens_id), 7 (PIP), 8 (TPM), 9 (distance), 
# 10-12 (enhancer loc), 13 (enhancer cell type), 14 (enhancer target gene hgnc), 15 (enhancer score)
rule intersect_variants_predictions:
	input:
		predictionsSorted = os.path.join(config["outDir"], "{method}", "{biosample}", "enhancerPredictions.sorted.bed.gz"),
		filteredGTExVariantsFinal = os.path.join(config["outDir"], "{method}", "GTExVariants.filteredForMethod.tsv")
	params: 
		outDir = config["outDir"],
		chrSizes = config["chrSizes"]
	output:
		variantsPredictionsInt = os.path.join(config["outDir"], "{method}", "{biosample}", "GTExVariants-enhancerPredictionsInt.tsv.gz")
	conda: 
		os.path.join(config["envDir"], "eQTLEnv.yml")
	shell:
		"""
		set +o pipefail;
		
		zcat {input.predictionsSorted} | bedtools intersect -wa -wb -sorted -a {input.filteredGTExVariantsFinal} -b stdin -g {params.chrSizes} | gzip > {output.variantsPredictionsInt}
		"""

# output columns: 1-3 (loc), 4 (rsID), 5-7 (enhancer loc), 8 (enhancer cell type), 9 (enhancer target gene hgnc), 10 (enhancer score)
rule intersect_bg_variants_predictions:
	input:
		predictionsSorted = os.path.join(config["outDir"], "{method}", "{biosample}", "enhancerPredictions.sorted.bed.gz"),
		commonVarDistalNoncoding = os.path.join(config["outDir"], "generalReference", "distalNoncodingBackgroundSNPs.bed.gz")
	params:
		outDir = config["outDir"],
		chrSizes = config["chrSizes"]
	output:
		commonVarPredictionsInt = os.path.join(config["outDir"], "{method}", "{biosample}", "distalNoncodingBackgroundSNPs-enhancerPredictionsInt.tsv.gz")
	conda: 
		os.path.join(config["envDir"], "eQTLEnv.yml")
	shell:
		"""
		set +o pipefail;
		
		zcat {input.predictionsSorted} | bedtools intersect -wa -wb -sorted -a {input.commonVarDistalNoncoding} -b stdin -g {params.chrSizes} | gzip > {output.commonVarPredictionsInt}
		"""

# generate threshold span based on quantiles of interescting predictions
rule generate_quantile_threshold_span:
	input:
		# read in internally
	params:
		methods_config = config["methodsTable"],
		nSteps = config["nThresholdSteps"],
		outDir = config["outDir"],
	output:
		outFile = os.path.join(config["outDir"], "{method}", "thresholdSpan.tsv")
	conda:
		os.path.join(config["envDir"], "eQTLEnv.yml")	
	script:
		os.path.join(config["codeDir"], "generate_quantile_threshold_span.R")

