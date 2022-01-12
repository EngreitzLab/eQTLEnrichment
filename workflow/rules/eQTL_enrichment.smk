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
			Rscript {params.codeDir}/filter_to_ABC_genes.R --input {input.GTExGeneUniverse} --col 4 --genes {input.methodGeneUniverse} --id hgnc > {output.geneUniverse}
			"""

# return list of biosamples for method from biosample key
rule get_biosamples:
	input: 
		biosampleKey = lambda wildcards: methods_config.loc[wildcards.method, "sampleKey"]
	output:
		samples = os.path.join(config["outDir"], "{method}", "biosampleList.tsv")
	conda: 
		os.path.join(config["envDir"], "eQTLEnv.yml")
	script: 
		os.path.join(config["codeDir"], "get_biosamples.R")
		
# sort enhancer predictions by chromosome & start location & filter to gene universe, get list of samples
rule sort_predictions:
	input:
		predFile = lambda wildcards: methods_config.loc[wildcards.method, "predFiles"][wildcards.biosample],
		geneUniverse = os.path.join(config["outDir"], "{method}", "geneUniverse.tsv")
	params:
		codeDir = config["codeDir"],
		outDir = config["outDir"],
		chrSizes = config["chrSizes"]
	output:
		predictionsSorted = os.path.join(config["outDir"], "{method}", "{biosample}", "enhancerPredictions.sorted.bed.gz"),
	conda: 
		os.path.join(config["envDir"], "eQTLEnv.yml")
	shell:
		"""
		set +o pipefail;
			
		# sort predictions file
		zcat {input.predFile} | csvtk cut -t -f chr,start,end,CellType,TargetGene,Score | sed 1d | bedtools sort -i stdin -faidx {params.chrSizes} > {params.outDir}/{wildcards.method}/{wildcards.biosample}/temp.sortedPred.tsv
			
		# filter predictions to gene universe
		Rscript {params.codeDir}/filter_to_ABC_genes.R --input {params.outDir}/{wildcards.method}/{wildcards.biosample}/temp.sortedPred.tsv --col 5 --genes {input.geneUniverse} --id hgnc | gzip > {output.predictionsSorted}
			
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
		chrSizes = config["chrSizes"]
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
			awk '$4=="ABC" || $4=="AllPeaks" || $4=="Other" || $4=="OtherIntron"' {input.partition} | bedtools sort -i stdin -faidx {params.chrSizes} > {output.partitionDistalNoncoding}
			
			# filter GTEx variants by SUSIE, credible set, PIP0.5; print set of columns: 1-3 (loc), 4 (hgID), 5 (tissue), 6 (ens_id), 7 (PIP); then filter to distal noncoding
			zcat {input.GTExVariants} | awk '$16>=0.5 && $17 != -1  && $9 == "SUSIE"' | cut -f1-4,10,11,16 | bedtools sort -i stdin -faidx {params.chrSizes} | uniq | bedtools intersect -wa -sorted -a stdin -b {output.partitionDistalNoncoding} -g {params.chrSizes} | gzip > {output.filteredGTExVar}

			# filter common variants to distal noncoding
			cat {input.commonVar} | bedtools intersect -wa -sorted -a stdin -b {output.partitionDistalNoncoding} -g {params.chrSizes}| gzip > {output.commonVarDistalNoncoding}
			"""

# filter variants by expression (based on what tissue they are in and GTEx expression data)
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

# filter variants to by gene universe for each method
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
	
# intersect predictions for each threshold with GTEx variants
rule intersect_variants_predictions:
	input:
		predictionsSorted = os.path.join(config["outDir"], "{method}", "{biosample}", "enhancerPredictions.thresholded.{threshold}.bed.gz"),
		filteredGTExVariantsFinal = os.path.join(config["outDir"], "{method}", "GTExVariants.filteredForMethod.tsv")
	params: 
		outDir = config["outDir"],
		chrSizes = config["chrSizes"]
	output:
		variantsPredictionsInt = os.path.join(config["outDir"], "{method}", "{biosample}", "GTExVariants-enhancerPredictionsInt.{threshold}.tsv.gz")
	conda: 
		os.path.join(config["envDir"], "eQTLEnv.yml")
	shell:
		"""
		set +o pipefail;
		
		# columns of output: 1-3 (variant loc), 4 (variant hgID), 5 (variant tissue), 6 (eGene), 7 (variant PIP), 8 (eGene TPM), 9-11 (enhancer loc), 12 (enhancer cell type), 13 (enhancer target gene hgnc), 14 (enhancer score)
		zcat {input.predictionsSorted} | bedtools intersect -wa -wb -sorted -a {input.filteredGTExVariantsFinal} -b stdin -g {params.chrSizes} | gzip > {output.variantsPredictionsInt}
		"""

# compute count matrix (number of eQTLSs in predicted enhancers at every GTEx tissue/biosample intersection)
# in second rule all	
rule compute_count_matrix:
	input:
		samples = os.path.join(config["outDir"], "{method}", "biosampleList.tsv")
	params:
		outDir = config["outDir"],
		GTExTissues=config["GTExTissues"],
	output: 
		countMatrix = os.path.join(config["outDir"], "{method}", "count_matrix.{threshold}.tsv")
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
			
# computer number of common variants overlapping enhancers in each biosample for each threshold
rule compute_common_var_overlap:
	input:
		commonVarDistalNoncoding = os.path.join(config["outDir"], "generalReference", "distalNoncoding.bg.SNPs.bed.gz"),
		predictionsSorted = os.path.join(config["outDir"], "{method}", "{biosample}", "enhancerPredictions.thresholded.{threshold}.bed.gz"),
		samples = os.path.join(config["outDir"], "{method}", "biosampleList.tsv")
	params:
		outDir = config["outDir"],
		codeDir = config["codeDir"],
		chrSizes = config["chrSizes"]
	output:
		commonVarPerBiosample = os.path.join(config["outDir"], "{method}", "{biosample}", "commonVarPerBiosample.{threshold}.tsv")
	conda: 
		os.path.join(config["envDir"], "eQTLEnv.yml")
	shell:
		"""			
		set +o pipefail;
					
		zcat {input.predictionsSorted} | bedtools intersect -wa -wb -sorted -g {params.chrSizes} -a stdin -b {input.commonVarDistalNoncoding}| cut -f 4 | sort | uniq -c > {output.commonVarPerBiosample}
		"""

# generate matrix with enrichment values for each GTEx tissue/biosample intersection for each threshold
rule compute_enrichment_matrix:
	input: 
		countMatrix = os.path.join(config["outDir"], "{method}", "count_matrix.{threshold}.tsv"),
		variantsPerGTExTissue = os.path.join(config["outDir"], "{method}", "nVariantsPerGTExTissue.tsv"),
		#  read in within method, see note in count matrix
		samples = os.path.join(config["outDir"],"{method}", "biosampleList.tsv"),
	params:
		codeDir = config["codeDir"],
		outDir = config["outDir"],
		GTExTissues = config["GTExTissues"],
		sampleKey = lambda wildcards: methods_config.loc[wildcards.method, "sampleKey"],
	output: 
		enrichmentTable = os.path.join(config["outDir"], "{method}", "enrichmentTable.{threshold}.tsv")
	conda: 
		os.path.join(config["envDir"], "eQTLEnv.yml")
	script: 
		os.path.join(config["codeDir"], "counts_to_enrichment.R")
