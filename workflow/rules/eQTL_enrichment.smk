
# get list of genes considered by GTEx and by the prediction method
# run once per method
rule make_gene_universes:
	input: 
 		methodGeneUniverse = lambda wildcards: config["geneUniverse"][wildcards.method],
		GTExGeneUniverse = config["GTExGeneUniverse"]
	params:
		outDir = config["outDir"],
		codeDir = config["codeDir"]
	output:
		geneUniverse = os.path.join(config["outDir"], "{method}", "geneUniverse.tsv")
	run:
		shell(
			"""
			set +o pipefail;
			
			# filter GTEx genes to method genes
			Rscript {params.codeDir}/filter_to_ABC_genes.R --input {input.GTExGeneUniverse} --col 4 --genes {input.methodGeneUniverse} --id hgnc > {output.geneUniverse}
			""")

# sort enhancer predictions by chromosome & start location & filter to gene universe, get list of samples
# run once per prediction method
rule sort_predictions:
	input:
		predFile = lambda wildcards: config["predFile"][wildcards.method],
		geneUniverse = os.path.join(config["outDir"], "{method}", "geneUniverse.tsv")
	params:
		codeDir = config["codeDir"],
		outDir = config["outDir"]
	output:
		predictionsSorted = os.path.join(config["outDir"], "{method}", "enhancerPredictions.sorted.bed.gz"),
		samples = os.path.join(config["outDir"], "{method}", "biosampleList.tsv")
	run:
		shell(
			"""
			set +o pipefail;
			
			# get list of samples alphabetically (in a column)
			zcat {input.predFile} | csvtk cut -t -f CellType | sed 1d | sort -k1,1 | uniq > {output.samples}
			
			# sort predictions file
			zcat {input.predFile} | csvtk cut -t -f chr,start,end,CellType,TargetGene | sed 1d | sort -k1,1 -k2,2n > {params.outDir}/{wildcards.method}/temp.sortedPred.tsv
			
			# filter predictions to gene universe
			Rscript {params.codeDir}/filter_to_ABC_genes.R --input {params.outDir}/{wildcards.method}/temp.sortedPred.tsv --col 4 --genes {input.geneUniverse} --id hgnc | gzip > {output.predictionsSorted}
			
			rm {params.outDir}/{wildcards.method}/temp.sortedPred.tsv
			
			""")

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
		commonVarDistalNoncoding = os.path.join(config["outDir"], "generalReference", "distalNoncoding.bg.SNPs.bed.gz"),	
	run:
		shell(
			"""
			set +o pipefail;
			
			# filter partition to distal noncoding
			awk '$4=="ABC" || $4=="AllPeaks" || $4=="Other" || $4=="OtherIntron"' {input.partition} | sort -k1,1 -k2,2n > {output.partitionDistalNoncoding}
			
			# filter GTEx variants by SUSIE, credible set, PIP0.5; print set of columns: 1-3 (loc), 4 (hgID), 5 (tissue), 6 (ens_id), 7 (PIP); then filter to distal noncoding
			zcat {input.GTExVariants} | awk '$16>=0.5 && $17 != -1  && $9 == "SUSIE"' | cut -f1-4,10,11,16 | sort -k1,1 -k2,2n | uniq | bedtools intersect -wa -sorted -a stdin -b {output.partitionDistalNoncoding} | gzip > {output.filteredGTExVar}

			# filter common variants to distal noncoding
			cat {input.commonVar} | bedtools intersect -wa -sorted -a stdin -b {output.partitionDistalNoncoding} | gzip > {output.commonVarDistalNoncoding}
			""")

# filter variants by expression (once overall)
rule filter_variants_by_expression:
	input: 
		exprData = config["GTExExpression"],
		filteredGTExVar = os.path.join(config["outDir"], "generalReference", "GTExVariants.filtered.PIP0.5.distalNoncoding.tsv.gz"),
	params:
		outDir = config["outDir"]
	output:
		GTExExpressedGenes = os.path.join(config["outDir"], "generalReference", "GTExVariants.filtered.PIP0.5.distalNoncoding.expressed.tsv")
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
	run:
		shell(
			"""
			## filter this method's gene universe
			Rscript {params.codeDir}/filter_to_ABC_genes.R --input {input.GTExExpressedGenes} --col 6 --genes {input.geneUniverse} --id hgnc > {output.filteredGTExVariantsFinal}
			
			""")	
# intersect predictions and GTEx variants
# run once per  method
rule intersect_variants_predictions:
	input:
		predictionsSorted = os.path.join(config["outDir"], "{method}", "enhancerPredictions.sorted.bed.gz"),
		filteredGTExVariantsFinal = os.path.join(config["outDir"], "{method}", "GTExVariants.filteredForMethod.tsv")
	params: 
		outDir = config["outDir"]
	output:
		variantsPredictionsInt = os.path.join(config["outDir"], "{method}", "GTExVariants-enhancerPredictionsInt.tsv.gz")
	run:
		shell(
			"""
			# columns of output: 1-3 (variant loc), 4 (variant hgID), 5 (variant tissue), 6 (egene), 7 (variant PIP), 8 (eGene TPM), 9-11 (enhancer loc), 12 (enhancer cell type), 13 (enhancer target gene hgnc)
			zcat {input.predictionsSorted} | bedtools intersect -wa -wb -sorted -a {input.filteredGTExVariantsFinal} -b stdin | gzip > {output.variantsPredictionsInt}
			""")
        
# make count matrix
# once per method
rule compute_count_matrix_old:
	input:
		variantsPredictionsInt = os.path.join(config["outDir"], "{method}", "GTExVariants-enhancerPredictionsInt.tsv.gz"),
		samples = os.path.join(config["outDir"], "{method}", "biosampleList.tsv")
	params:
		outDir = config["outDir"],
		GTExTissues=config["GTExTissues"]
	output: 
		countMatrix = os.path.join(config["outDir"], "{method}", "countMatrixOld.tsv")
	run:
		shell(
			"""			
			set +o pipefail;
			
			biosamples=$(awk 'BEGIN {{ ORS=" " }} {{ print }}' {input.samples})
			for sample in $biosamples
			do
			# columns of output: 1-3 (variant loc), 4 (variant hgID), 5 (variant tissue), 6 (egene), 7 (variant PIP), 8 (eGene TPM), 9-11 (enhancer loc), 12 (enhancer cell type), 13 (enhancer target gene hgnc)
				# filter for cell type, select just variant hgID and tissue, deduplicate (through uniq)
				# print counts per GTEx tissue, sort (alphabetically by first tissue), print only count, then print all numbers in a row in matrix
				zcat {input.variantsPredictionsInt} | awk -v awkvar=$sample '$12==awkvar' | cut -f 4,5 | sort -k1 | uniq | awk '{{a[$2]++}} END {{for(i in a) print i" "a[i]}}' | sort -k1 | cut -d' ' -f2 | awk 'BEGIN {{ ORS = "\\t" }} {{ print }}' >> {output.countMatrix}
				printf "\\n" >> {output.countMatrix}
			done
			""")
			
rule compute_count_matrix:
	input:
		variantsPredictionsInt = os.path.join(config["outDir"], "{method}", "GTExVariants-enhancerPredictionsInt.tsv.gz"),
		samples = os.path.join(config["outDir"], "{method}", "biosampleList.tsv")
	params:
		outDir = config["outDir"],
		GTExTissues=config["GTExTissues"]
	output: 
		countMatrix = os.path.join(config["outDir"], "{method}", "countMatrix.tsv")
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
	run:
		shell(
			"""
			set +o pipefail;
			# get counts per GTEx tissue, columns of GTEx var: 1-3 (loc), 4 (hgID), 5 (tissue), 6 (gene ens id), 7 (PIP)

			for tissue in {params.GTExTissues}
			do
				printf $tissue"\\t" >> {output.variantsPerGTExTissue}
				cat {input.filteredGTExVariantsFinal} | awk -v awkvar=$tissue '$5==awkvar' | wc -l >> {output.variantsPerGTExTissue}
			done
			""")
			
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
		commonVarPerBiosample = os.path.join(config["outDir"], "{method}", "commonVarPerBiosample.tsv"),
		commonVarPredictionsInt = os.path.join(config["outDir"], "{method}", "commonVarPredictionsInt.tsv")
	run:
		shell(
			"""			
			set +o pipefail;
			
			# intersect common variants with predictions
			# output columns: (1-3) common variant loc, (4) common variant rsid, (5-7) enhancer loc, (8) cell type, (9) target gene
			zcat {input.commonVarDistalNoncoding} > {params.outDir}/{wildcards.method}/tempCommonVar.tsv
			zcat {input.predictionsSorted} | bedtools intersect -wa -wb -sorted -a {params.outDir}/{wildcards.method}/tempCommonVar.tsv -b stdin > {output.commonVarPredictionsInt}
			
			# get common variants: biosample/celltype followed by count
			biosamples=$(awk 'BEGIN {{ ORS=" " }} {{ print }}' {input.samples})
			for sample in $biosamples
			do	
				printf $sample'\\t' >> {output.commonVarPerBiosample}
				cat {output.commonVarPredictionsInt} | awk -v awkvar=$sample '$8==awkvar' | wc -l >> {output.commonVarPerBiosample}
			done
			
			rm {params.outDir}/{wildcards.method}/tempCommonVar.tsv
			""")
		
# compute number of base pairs in each enhancer set
rule compute_enhancer_set_size:
	input:
		predictionsSorted = os.path.join(config["outDir"], "{method}", "enhancerPredictions.sorted.bed.gz"),
		samples = os.path.join(config["outDir"],"{method}", "biosampleList.tsv")
	params:
	output:
		basesPerEnhancerSet = os.path.join(config["outDir"], "{method}", "basesPerEnhancerSet.tsv")
	run:
		shell(
			"""			
			set +o pipefail;
			
			biosamples=$(awk 'BEGIN {{ ORS=" " }} {{ print }}' {input.samples})
			
			for sample in $biosamples
			do
				# make list of sizes of enhancer sets for each biosample
				printf $sample"\\t" >> {output.basesPerEnhancerSet}
				zcat {input.predictionsSorted} | awk -v awkvar=$sample '$4==awkvar' | sort -k1,1 -k2,2n | bedtools merge -i stdin | awk 'BEGIN {{FS=OFS="\\t"}} {{print $3-$2}}' | awk '{{s+=$1}}END{{print s}}' >> {output.basesPerEnhancerSet}
			done
			""")

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
		sampleKey = lambda wildcards: config["enrichment_heatmaps"]["sampleKey"][wildcards.method],
		sampleID = lambda wildcards: config["enrichment_heatmaps"]["sampleID"][wildcards.method],
		sampleName = lambda wildcards: config["enrichment_heatmaps"]["sampleName"][wildcards.method]
	output: 
		enrichmentTable = os.path.join(config["outDir"], "{method}", "enrichmentTable.tsv")
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
		sampleKey = lambda wildcards: config["enrichment_heatmaps"]["sampleKey"][wildcards.method],
		sampleID = lambda wildcards: config["enrichment_heatmaps"]["sampleID"][wildcards.method],
		sampleName = lambda wildcards: config["enrichment_heatmaps"]["sampleName"][wildcards.method],
		sampleCategory = lambda wildcards: config["enrichment_heatmaps"]["sampleCategory"][wildcards.method]
	output: 
		heatmapFull = os.path.join(config["outDir"], "{method}", "enrichmentHeatmap.full.pdf"),
		heatmapAggregated = os.path.join(config["outDir"], "{method}", "enrichmentHeatmap.aggregated.pdf")
	run:
		shell(
			"""	
			set +o pipefail;
			
			# plot full heat map
			Rscript {params.codeDir}/plot_enrichment_heatmap.R --table {input.enrichmentTable} --outfile {output.heatmapFull} --samplekey {params.sampleKey} --cellid {params.sampleID} --category {params.sampleName} --enhancersizes {params.basesPerEnhancerSet}
			
			# plot aggregated heat map
			Rscript {params.codeDir}/plot_enrichment_heatmap.R --table {input.enrichmentTable} --outfile {output.heatmapAggregated} --samplekey {params.sampleKey} --cellid {params.sampleID} --category {params.sampleCategory} --enhancersizes {params.basesPerEnhancerSet}
			
			""")
			
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
	script: 
		os.path.join(config["codeDir"], "enrichment_report.Rmd")