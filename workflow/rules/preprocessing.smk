# intersect genes considered by variants and prediction method
rule make_gene_universes:
	input: 
 		methodGeneUniverse = lambda wildcards: methods_config.loc[wildcards.method, "geneUniverse"],
		GTExGeneUniverse = config["GTExGeneUniverse"]
	output:
		geneUniverse = temp(os.path.join(config["outDir"], "{method}", "intermediate", "geneUniverse.bed.gz"))
	conda: 
		os.path.join(config["envDir"], "eQTLEnv.yml")
	resources:
		mem_mb = determine_mem_mb
	shell:
			"""
			set +o pipefail;
			
			# return variant universe genes that are also in method gene universe
			awk 'NR==FNR{{names[$4]; next}} $4 in names' {input.GTExGeneUniverse} {input.methodGeneUniverse} | gzip > {output.geneUniverse}

			"""

# sort enhancer predictions by chromosome & start location & filter to gene universe, invert scores if necessary
# return file with (1-3) loc, (4) biosample, (5) TargetGene, (6) score (no header)
rule process_predictions:
	input:
		predFile = lambda wildcards: methods_config.loc[wildcards.method, "predFiles"][wildcards.biosample],
		geneUniverse = os.path.join(config["outDir"], "{method}", "intermediate", "geneUniverse.bed.gz")
	params:
		codeDir = config["codeDir"],
		outDir = config["outDir"],
		chrSizes = config["chrSizes"],
		scoreCol = lambda wildcards: methods_config.loc[wildcards.method, "score_col"],
		inversePred = lambda wildcards: methods_config.loc[wildcards.method, "inverse_predictor"]
	output:
		predictionsSorted = temp(os.path.join(config["outDir"], "{method}", "biosamples", "{biosample}", "enhancerPredictions.sorted.bed.gz"))
	resources:
		mem_mb = determine_mem_mb
	conda: 
		os.path.join(config["envDir"], "eQTLEnv.yml")
	shell:
		"""
		set +o pipefail;
			
		# sort predictions file: remove # from header,select columns,remove header, remove rows with blanks
        if [[ {input.predFile} == *.gz ]]
        then
		    zcat {input.predFile} | awk 'NR==1{{sub(/^#*/, "")}}1' | csvtk cut -t -f chr,start,end,TargetGene,{params.scoreCol} | sed 1d | awk 'NF==5{{print}}{{}}' | bedtools sort -i stdin -faidx {params.chrSizes} > {params.outDir}/{wildcards.method}/biosamples/{wildcards.biosample}/temp.sortedPred.tsv
		else
            cat {input.predFile} | awk 'NR==1{{sub(/^#*/, "")}}1' | csvtk cut -t -f chr,start,end,TargetGene,{params.scoreCol} | sed 1d | awk 'NF==5{{print}}{{}}' | bedtools sort -i stdin -faidx {params.chrSizes} > {params.outDir}/{wildcards.method}/biosamples/{wildcards.biosample}/temp.sortedPred.tsv
        fi

		# invert score if inverted predictor and filter to gene universe and set biosample column
		Rscript {params.codeDir}/preprocessing/process_predictions.R --input {params.outDir}/{wildcards.method}/biosamples/{wildcards.biosample}/temp.sortedPred.tsv  --genes {input.geneUniverse} --biosample {wildcards.biosample} --invert {params.inversePred}  | gzip > {output.predictionsSorted}

		rm {params.outDir}/{wildcards.method}/biosamples/{wildcards.biosample}/temp.sortedPred.tsv
			
		"""

# filter eQTLs by PIP and to distal noncoding; filter common variants to distal noncoding
rule filter_all_variants:
	input:
		GTExVariants = config["eQTLVariants"],
		commonVar = config["bgVariants"],
		partition = config["partition"]
	params: 
		chrSizes = config["chrSizes"],
		thresholdPIP = config["thresholdPIP"]
	output: 
		filteredGTExVar = temp(os.path.join(config["outDir"], "variants", "GTExVariants.PIPfilt.distalNoncoding.tsv.gz")),
		partitionDistalNoncoding = temp(os.path.join(config["outDir"], "variants", "Partition.distalNoncoding.bed")),
		commonVarDistalNoncoding = temp(os.path.join(config["outDir"], "variants", "distalNoncodingBackgroundSNPs.bed.gz"))
	resources:
		mem_mb = determine_mem_mb
	conda: 
		os.path.join(config["envDir"], "eQTLEnv.yml")	
	shell:
			"""
			set +o pipefail;
			
			# filter partition to distal noncoding
			awk '$4=="ABC" || $4=="AllPeaks" || $4=="Other" || $4=="OtherIntron"' {input.partition} | bedtools sort -i stdin -faidx {params.chrSizes} > {output.partitionDistalNoncoding}

			# variant file: reorder columns, filter by PIP, filter to distal noncoding
			zcat {input.GTExVariants} | csvtk cut -t -f chr,start,end,varID_hg38,gene_hgnc,tissue,pip | sed 1d |  awk '$7>={params.thresholdPIP}' | bedtools sort -i stdin -faidx {params.chrSizes} | uniq | bedtools intersect -wa -sorted -a stdin -b {output.partitionDistalNoncoding}  -g {params.chrSizes} | gzip > {output.filteredGTExVar}

			# filter common variants to distal noncoding
			cat {input.commonVar} | bedtools sort -i stdin -faidx {params.chrSizes} | bedtools intersect -wa -sorted -a stdin -b {output.partitionDistalNoncoding} -g {params.chrSizes} | gzip > {output.commonVarDistalNoncoding}
			"""

rule bg_variant_count:
	input:
		commonVarDistalNoncoding = os.path.join(config["outDir"], "variants", "distalNoncodingBackgroundSNPs.bed.gz"),
	output:
		commonVarCounts= os.path.join(config["outDir"], "variants", "distalNoncodingBackgroundSNPCount.txt")
	shell:
		"""
		set +o pipefail;

		zcat {input.commonVarDistalNoncoding} | wc -l > {output.commonVarCounts}

		"""

# add eVariant - eGene TSS variant distance to variant file
# columns: 1-3 (loc), 4 (variantID), 5 (gene), 6 (tissue), 7 (PIP),  8 (distance group)
rule add_distance_to_variants:
	input:
		filteredGTExVariants = os.path.join(config["outDir"], "variants", "GTExVariants.PIPfilt.distalNoncoding.tsv.gz"),
	params:
		TSS = config['TSS'],
		distances = config["distances"]
	output:
		GTExVariantsDistance = temp(os.path.join(config["outDir"], "variants", "GTExVariants.PIPfilt.distalNoncoding.withDistance.tsv.gz"))
	resources:
		mem_mb = determine_mem_mb
	conda:
		os.path.join(config["envDir"], "eQTLEnv.yml")
	script:
		os.path.join(config["codeDir"], "preprocessing", "add_distance_to_variants.R")


# filter variants to by gene universe for each method
# columns: 1-3 (loc), 4 (variantID), 5 (gene), 6 (tissue), 7 (PIP),  8 (distance group)
rule filter_variants_to_gene_universe:
	input:
		GTExVariantsDistance = os.path.join(config["outDir"], "variants", "GTExVariants.PIPfilt.distalNoncoding.withDistance.tsv.gz"),
		geneUniverse = os.path.join(config["outDir"], "{method}", "intermediate", "geneUniverse.bed.gz")
	output:
		filteredGTExVariantsFinal = temp(os.path.join(config["outDir"], "{method}", "intermediate", "GTExVariants.filteredForMethod.tsv.gz"))
	resources:
		mem_mb = determine_mem_mb
	shell:
		"""
		set +o pipefail;

		# filter variants based on gene universe
		awk 'NR==FNR{{names[$4]; next}} $5 in names' <(zcat {input.geneUniverse}) <(zcat {input.GTExVariantsDistance}) | gzip > {output.filteredGTExVariantsFinal}
			
		"""

# intersect predictions for each threshold with GTEx variants
# output columns:  columns: 1-3 (loc), 4 (variantID), 5 (gene), 6 (tissue), 7 (PIP),  8 (distance group)
# 9-11 (enhancer loc), 12 (enhancer cell type), 13 (enhancer target gene hgnc), 14 (enhancer score)
rule intersect_variants_predictions:
	input:
		predictionsSorted = os.path.join(config["outDir"], "{method}", "biosamples", "{biosample}", "enhancerPredictions.sorted.bed.gz"),
		filteredGTExVariantsFinal = os.path.join(config["outDir"], "{method}", "intermediate", "GTExVariants.filteredForMethod.tsv.gz")
	params: 
		chrSizes = config["chrSizes"]
	output:
		variantsPredictionsInt = temp(os.path.join(config["outDir"], "{method}", "biosamples", "{biosample}", "GTExVariants-enhancerPredictionsInt.tsv.gz"))
	resources:
		mem_mb = determine_mem_mb
	conda: 
		os.path.join(config["envDir"], "eQTLEnv.yml")
	shell:
		"""
		set +o pipefail;

		zcat {input.predictionsSorted} | bedtools intersect -wa -wb -sorted -a <(zcat {input.filteredGTExVariantsFinal}) -b stdin -g {params.chrSizes} | gzip > {output.variantsPredictionsInt}

		"""

# output columns: 1-3 (loc), 4 (rsID), 5-7 (enhancer loc), 8 (enhancer cell type), 9 (enhancer target gene hgnc), 10 (enhancer score)
rule intersect_bg_variants_predictions:
	input:
		predictionsSorted = os.path.join(config["outDir"], "{method}", "biosamples", "{biosample}", "enhancerPredictions.sorted.bed.gz"),
		commonVarDistalNoncoding = os.path.join(config["outDir"], "variants", "distalNoncodingBackgroundSNPs.bed.gz")
	params:
		chrSizes = config["chrSizes"]
	output:
		commonVarPredictionsInt = temp(os.path.join(config["outDir"], "{method}", "biosamples", "{biosample}", "distalNoncodingBackgroundSNPs-enhancerPredictionsInt.tsv.gz")),
	resources:
		mem_mb = determine_mem_mb
	conda: 
		os.path.join(config["envDir"], "eQTLEnv.yml")
	shell:
		"""
		set +o pipefail;
		
		zcat {input.predictionsSorted} | bedtools intersect -wa -wb -sorted -a <(zcat {input.commonVarDistalNoncoding}) -b stdin -g {params.chrSizes} | gzip > {output.commonVarPredictionsInt}
		"""

# make tsv with columns "tissue, biosample"
rule biosample_tissue_maps:
	params:
		tissues = lambda wildcards: methods_config.loc[wildcards.method, "GTExTissue_map"],
		biosamples = lambda wildcards: methods_config.loc[wildcards.method, "biosample_map"]
	output:
		map = os.path.join(config["outDir"], "{method}", "intermediate", "GTExTissueBiosampleMap.tsv")
	resources:
		mem_mb = determine_mem_mb
	conda: 
		os.path.join(config["envDir"], "eQTLEnv.yml")
	script:
		os.path.join(config["codeDir"], "preprocessing", "make_biosample_tissue_map.R")

# generate threshold span based on quantiles of interescting predictions
rule generate_quantile_threshold_span:
	input:
		map = os.path.join(config["outDir"], "{method}", "intermediate", "GTExTissueBiosampleMap.tsv"),
		varInt = lambda wildcards: [os.path.join(config["outDir"], wildcards.method, "biosamples", biosample, "GTExVariants-enhancerPredictionsInt.tsv.gz") for biosample in methods_config.loc[wildcards.method, "biosamples"]]
	params:
		nSteps = config["nThresholdSteps"],
		biosamples = lambda wildcards: methods_config.loc[wildcards.method, "biosamples"],
		binary = lambda wildcards: methods_config.loc[wildcards.method, "boolean"]
	output:
		outFile = temp(os.path.join(config["outDir"], "{method}", "intermediate", "thresholdSpan.tsv"))
	resources:
		mem_mb = determine_mem_mb
	conda:
		os.path.join(config["envDir"], "eQTLEnv.yml")	
	script:
		os.path.join(config["codeDir"], "preprocessing", "generate_quantile_threshold_span.R")

