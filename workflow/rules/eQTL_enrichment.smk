# load prediction method config file
methods_config_file = config["methodsTable"]
methods_config = pd.read_table(methods_config_file).set_index("method", drop=False)
methods_config["GTExTissue_map"] = methods_config["GTExTissue_map"].apply(eval)
methods_config["biosample_map"] = methods_config["biosample_map"].apply(eval)

# get list of genes considered by GTEx and by the prediction method
# run once per method
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
			Rscript {params.codeDir}/filter_to_ABC_genes.R --input {input.GTExGeneUniverse} --col 4 --genes {input.methodGeneUniverse} --id hgnc > {output.geneUniverse}
			"""

# sort enhancer predictions by chromosome & start location & filter to gene universe, get list of samples
# run once per prediction method
rule sort_predictions:
	input:
		predFile = lambda wildcards: methods_config.loc[wildcards.method, "predFile"],
		geneUniverse = os.path.join(config["outDir"], "{method}", "geneUniverse.tsv")
	params:
		codeDir = config["codeDir"],
		outDir = config["outDir"]
	output:
		predictionsSorted = os.path.join(config["outDir"], "{method}", "enhancerPredictions.sorted.bed.gz"),
		samples = os.path.join(config["outDir"], "{method}", "biosampleList.tsv")
	conda: 
		os.path.join(config["envDir"], "eQTLEnv.yml")
	shell:
		"""
		set +o pipefail;
			
		# get list of samples alphabetically (in a column)
		zcat {input.predFile} | csvtk cut -t -f CellType | sed 1d | sort -k1,1 | uniq > {output.samples}
			
		# sort predictions file
		zcat {input.predFile} | csvtk cut -t -f chr,start,end,CellType,TargetGene | sed 1d | sort -k1,1 -k2,2n > {params.outDir}/{wildcards.method}/temp.sortedPred.tsv
			
		# filter predictions to gene universe
		Rscript {params.codeDir}/filter_to_ABC_genes.R --input {params.outDir}/{wildcards.method}/temp.sortedPred.tsv --col 5 --genes {input.geneUniverse} --id hgnc | gzip > {output.predictionsSorted}
			
		rm {params.outDir}/{wildcards.method}/temp.sortedPred.tsv
			
		"""

# filter GTEx variants by PIP, credible set, and to distal noncoding genes; convert ensembl id to hgnc
# run once overall
rule filter_all_variants:
	input:
		GTExVariants = config["GTExVariants"],
		commonVar = config["bgVariants"],
		partition = config["partition"]
	params: 
		codeDir = config["codeDir"],
		outDir = config["outDir"]
	output: 
		filteredGTExVar = os.path.join(config["outDir"], "generalReference", "GTExVariants.filtered.PIP0.5.distalNoncoding.tsv.gz"),
		partitionDistalNoncoding = os.path.join(config["outDir"], "generalReference", "Partition.distalNoncoding.bed"),
		commonVarDistalNoncoding = os.path.join(config["outDir"], "generalReference", "distalNoncoding.bg.SNPs.bed.gz")
	conda: 
		os.path.join(config["envDir"], "eQTLEnv.yml")	
	shell:
			"""
			set +o pipefail;
			
			# filter partition to distal noncoding
			awk '$4=="ABC" || $4=="AllPeaks" || $4=="Other" || $4=="OtherIntron"' {input.partition} | sort -k1,1 -k2,2n > {output.partitionDistalNoncoding}
			
			# filter GTEx variants by SUSIE, credible set, PIP0.5; print set of columns: 1-3 (loc), 4 (hgID), 5 (tissue), 6 (ens_id), 7 (PIP); then filter to distal noncoding
			zcat {input.GTExVariants} | awk '$16>=0.5 && $17 != -1  && $9 == "SUSIE"' | cut -f1-4,10,11,16 | sort -k1,1 -k2,2n | uniq | bedtools intersect -wa -sorted -a stdin -b {output.partitionDistalNoncoding} | gzip > {output.filteredGTExVar}

			# filter common variants to distal noncoding
			cat {input.commonVar} | bedtools intersect -wa -sorted -a stdin -b {output.partitionDistalNoncoding} | gzip > {output.commonVarDistalNoncoding}
			"""

# filter variants by expression (once overall)
rule filter_variants_by_expression:
	input: 
		exprData = config["GTExExpression"],
		filteredGTExVar = os.path.join(config["outDir"], "generalReference", "GTExVariants.filtered.PIP0.5.distalNoncoding.tsv.gz"),
	params:
		outDir = config["outDir"]
	output:
		GTExExpressedGenes = os.path.join(config["outDir"], "generalReference", "GTExVariants.filtered.PIP0.5.distalNoncoding.expressed.tsv")
	conda: 
		os.path.join(config["envDir"], "eQTLEnv.yml")
	script: 
		os.path.join(config["codeDir"], "GTEx_expression.R")

# filter to by gene universe (per method)
rule filter_variants_to_gene_universe:
	input:
		GTExExpressedGenes = os.path.join(config["outDir"], "generalReference", "GTExVariants.filtered.PIP0.5.distalNoncoding.expressed.tsv"),
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
			Rscript {params.codeDir}/filter_to_ABC_genes.R --input {input.GTExExpressedGenes} --col 6 --genes {input.geneUniverse} --id hgnc > {output.filteredGTExVariantsFinal}
			
		"""
	
# intersect predictions and GTEx variants
# run once per method
rule intersect_variants_predictions:
	input:
		predictionsSorted = os.path.join(config["outDir"], "{method}", "enhancerPredictions.sorted.bed.gz"),
		filteredGTExVariantsFinal = os.path.join(config["outDir"], "{method}", "GTExVariants.filteredForMethod.tsv")
	params: 
		outDir = config["outDir"]
	output:
		variantsPredictionsInt = os.path.join(config["outDir"], "{method}", "GTExVariants-enhancerPredictionsInt.tsv.gz")
	conda: 
		os.path.join(config["envDir"], "eQTLEnv.yml")
	shell:
		"""
		# columns of output: 1-3 (variant loc), 4 (variant hgID), 5 (variant tissue), 6 (eGene), 7 (variant PIP), 8 (eGene TPM), 9-11 (enhancer loc), 12 (enhancer cell type), 13 (enhancer target gene hgnc)
		zcat {input.predictionsSorted} | bedtools intersect -wa -wb -sorted -a {input.filteredGTExVariantsFinal} -b stdin | gzip > {output.variantsPredictionsInt}
		"""

# compute count matrix (number of eQTLSs in predicted enhancers at every GTEx tissue/biosample intersection)	
rule compute_count_matrix:
	input:
		variantsPredictionsInt = os.path.join(config["outDir"], "{method}", "GTExVariants-enhancerPredictionsInt.tsv.gz"),
		samples = os.path.join(config["outDir"], "{method}", "biosampleList.tsv")
	params:
		outDir = config["outDir"],
		GTExTissues=config["GTExTissues"]
	output: 
		countMatrix = os.path.join(config["outDir"], "{method}", "countMatrix.tsv")
	conda: 
		os.path.join(config["envDir"], "eQTLEnv.yml")
	script:
		os.path.join(config["codeDir"], "count_matrix.R")
	

# get number of variants per GTExTissue linked to expresed genes
rule get_variants_per_GTEx_tissue:
	input:
		filteredGTExVariantsFinal = os.path.join(config["outDir"], "{method}", "GTExVariants.filteredForMethod.tsv")
	params:
		outDir = config["outDir"],
		codeDir = config["codeDir"],
		GTExTissues = config["GTExTissues"]
	output: 
		variantsPerGTExTissue = os.path.join(config["outDir"], "{method}", "nVariantsPerGTExTissue.tsv")
	conda: 
		os.path.join(config["envDir"], "eQTLEnv.yml")
	shell:
		"""
		set +o pipefail;
		# get counts per GTEx tissue, columns of GTEx var: 1-3 (loc), 4 (hgID), 5 (tissue), 6 (gene ens id), 7 (PIP)

		for tissue in {params.GTExTissues}
		do
			printf $tissue"\\t" >> {output.variantsPerGTExTissue}
			cat {input.filteredGTExVariantsFinal} | awk -v awkvar=$tissue '$5==awkvar' | wc -l >> {output.variantsPerGTExTissue}
		done
		"""
			
# computer number of common variants overlapping enhancers in each biosample
# run once per method
rule compute_common_var_overlap:
	input:
		commonVarDistalNoncoding = os.path.join(config["outDir"], "generalReference", "distalNoncoding.bg.SNPs.bed.gz"),
		predictionsSorted = os.path.join(config["outDir"], "{method}", "enhancerPredictions.sorted.bed.gz"),
		samples = os.path.join(config["outDir"], "{method}", "biosampleList.tsv")
		
	params:
		outDir = config["outDir"],
		codeDir = config["codeDir"]
	output:
		commonVarPerBiosample = os.path.join(config["outDir"], "{method}", "commonVarPerBiosample.tsv")
	conda: 
		os.path.join(config["envDir"], "eQTLEnv.yml")
	shell:
		"""			
		set +o pipefail;
					
		zcat {input.predictionsSorted} | bedtools intersect -wa -wb -sorted -a stdin -b {input.commonVarDistalNoncoding}| cut -f 4 | sort | uniq -c > {output.commonVarPerBiosample}
		"""
		
# compute number of base pairs in each enhancer set
rule compute_enhancer_set_size:
	input:
		predictionsSorted = os.path.join(config["outDir"], "{method}", "enhancerPredictions.sorted.bed.gz"),
		samples = os.path.join(config["outDir"],"{method}", "biosampleList.tsv")
	params:
	output:
		basesPerEnhancerSet = os.path.join(config["outDir"], "{method}", "basesPerEnhancerSet.tsv")
	conda: 
		os.path.join(config["envDir"], "eQTLEnv.yml")
	shell:
		"""			
		set +o pipefail;
		
		biosamples=$(awk 'BEGIN {{ ORS=" " }} {{ print }}' {input.samples})
			
		for sample in $biosamples
		do
			# make list of sizes of enhancer sets for each biosample
			printf $sample"\\t" >> {output.basesPerEnhancerSet}
			zcat {input.predictionsSorted} | awk -v awkvar=$sample '$4==awkvar' | sort -k1,1 -k2,2n | bedtools merge -i stdin | awk 'BEGIN {{FS=OFS="\\t"}} {{print $3-$2}}' | awk '{{s+=$1}}END{{print s}}' >> {output.basesPerEnhancerSet}
		done
		"""

# generate matrix with enrichment values for each GTEx tissue/biosample intersection
# run once per method
rule compute_enrichment_matrix:
	input: 
		countMatrix = os.path.join(config["outDir"], "{method}","countMatrix.tsv"),
		variantsPerGTExTissue = os.path.join(config["outDir"], "{method}", "nVariantsPerGTExTissue.tsv"),
		commonVarPerBiosample = os.path.join(config["outDir"], "{method}", "commonVarPerBiosample.tsv"),
		samples = os.path.join(config["outDir"],"{method}", "biosampleList.tsv"),
	params:
		codeDir = config["codeDir"],
		outDir = config["outDir"],
		GTExTissues = config["GTExTissues"],
		sampleKey = lambda wildcards: methods_config.loc[wildcards.method, "sampleKey"],
	output: 
		enrichmentTable = os.path.join(config["outDir"], "{method}", "enrichmentTable.tsv")
	conda: 
		os.path.join(config["envDir"], "eQTLEnv.yml")
	script: 
		os.path.join(config["codeDir"], "counts_to_enrichment.R")

# plot aggregate and full heat maps of enrichment values
# run once per method
rule plot_enrichment_heatmaps:
	input: 
		enrichmentTable = os.path.join(config["outDir"], "{method}", "enrichmentTable.tsv")
	params:
		basesPerEnhancerSet = os.path.join(config["outDir"], "{method}", "basesPerEnhancerSet.tsv"),
		codeDir = config["codeDir"],
		sampleKey = lambda wildcards: methods_config.loc[wildcards.method, "sampleKey"]
	output: 
		heatmapFull = os.path.join(config["outDir"], "{method}", "enrichmentHeatmap.full.pdf"),
		heatmapAggregated = os.path.join(config["outDir"], "{method}", "enrichmentHeatmap.aggregated.pdf")
	conda: 
		os.path.join(config["envDir"], "eQTLEnv.yml")
	shell:
		"""	
		set +o pipefail;
			
		# plot full heat map
		Rscript {params.codeDir}/plot_enrichment_heatmap.R --table {input.enrichmentTable} --outfile {output.heatmapFull} --samplekey {params.sampleKey} --enhancersizes {params.basesPerEnhancerSet}  --useCategory "False"
			
		# plot aggregated heat map
		Rscript {params.codeDir}/plot_enrichment_heatmap.R --table {input.enrichmentTable} --outfile {output.heatmapAggregated} --samplekey {params.sampleKey} --enhancersizes {params.basesPerEnhancerSet} --useCategory "True"
			
		"""
			
# plot cdf and density graphs comparing enrichments across methods
# run once overall
rule plot_comparisons:
	input: 
		enrichmentTables = expand(os.path.join(config["outDir"], "{method}", "enrichmentTable.tsv"), method=config["methods"]),
	params:
		methods = config["methods"],
		codeDir = config["codeDir"],
		outDir = config["outDir"]
	output: 
		cdf = os.path.join(config["outDir"], "cdf.pdf"),
		density = os.path.join(config["outDir"], "density.pdf"),
		boxplot = os.path.join(config["outDir"], "boxplot.pdf")
	conda: 
		os.path.join(config["envDir"], "eQTLEnv.yml")
	script:
		os.path.join(config["codeDir"], "plot_comparisons.R")
			
# html report, run once overall
rule generate_report:
	input: 
		enrichmentTables = expand(os.path.join(config["outDir"],"{method}", "enrichmentTable.tsv"), method=config["methods"]),
		variantsPerGTExTissue = expand(os.path.join(config["outDir"], "{method}", "nVariantsPerGTExTissue.tsv"), method=config["methods"]),
		predictionMetrics = os.path.join(config["outDir"],"predictionMetrics.tsv"),
		sensitivitiesTable = os.path.join(config["outDir"],"sensitivitiesTable.tsv")
	params:
		names = config["methods"]
	output:
		enrichmentReport = os.path.join(config["outDir"], "enrichment_report.html")
	conda: 
		os.path.join(config["envDir"], "eQTLEnv.yml")
	script: 
		os.path.join(config["codeDir"], "enrichment_report.Rmd")