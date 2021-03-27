## rules for eQTL eGene prediction analysis

# filter GTEx variants by individual tissue
rule filter_variants_for_predictions:
	input: 
		GTExVariants = config["GTExVariants"],
		partitionDistalNoncoding = os.path.join(config["outDir"], "variantFilesForEnrichment/Partition.distalNoncoding.bed")
	params:
		codeDir = config["codeDir"],
		outDir = config["outDir"],
	output:
		filteredVariants = os.path.join(config["outDir"], "variantFilesForPrediction", "{GTExTissue}.filteredVariants.tsv.gz")
	run:
		shell("""
		set +o pipefail;
		
		# continue to use 0.5 PIP threshold?
		# just consider distal noncoding?
		# filter out non ABC genes?
		
		# column key: 1-3 = location, 4 = unique ID in hg19, 5 = hgnc ID, 6 = PIP, 7 = effect size, 8 = sd of effect size 
		zcat {input.GTExVariants} | awk -v var={wildcards.GTExTissue} '$16>=0.5 && $17 != -1  && $9 == "SUSIE" && $10==var' | cut -f 1,2,3,4,11,16,18,19 | bedtools intersect -wa -a stdin -b {input.partitionDistalNoncoding} | Rscript {params.codeDir}/get_hgnc_symbols.R --col 5 | sort -k1,1 -k2,2n | gzip > {output.filteredVariants}
		""")
		
# add columns for closest gene body, closest gene TSS, and whether eGene has its TSS within 100 kb for each variant
rule add_proximal_genes:
	input:
		filteredVariants = os.path.join(config["outDir"], "variantFilesForPrediction", "{GTExTissue}.filteredVariants.tsv.gz")
	params:
		chrSizes = config["chrSizes"],
		TSS = config["TSS"],
		genes = config["genes"],
		codeDir = config["codeDir"],
		outDir = config["outDir"]
	output:
		filteredVariantsProximal = os.path.join(config["outDir"], "variantFilesForPrediction", "{GTExTissue}.filteredVariants.proximalGenes.tsv.gz")
	run:
		shell("""
		set +o pipefail;
				
		# get list of distinct variants, expand variant +/- 100 kb, intersect with location of all TSS; end up with following columns: 1) variant unique ID, 2) gene with TSS within 100 kb
		zcat {input.filteredVariants} | cut -f 1-4 | sort -k 1,1 -k2,2n | uniq | bedtools slop -i stdin -b 100000 -g {params.chrSizes} | bedtools intersect -wa -wb -a stdin -b {params.TSS} | cut -f 4,8 > {params.outDir}/{wildcards.GTExTissue}.temp.100kbInterval.tsv

		# closest gene body (9), TSS (10), whether eGene has TSS within 100 kb (11)
		# how to handle cases where two genes or TSSs are equidistant? currently reporting only first (can do all, but leads to variants being counted multiple times)
		zcat {input.filteredVariants} | bedtools closest -t first -a stdin -b {params.genes} | cut -f 1-8,12 | bedtools closest -t first -a stdin -b {params.TSS} | cut -f  1-9,13 > {params.outDir}/{wildcards.GTExTissue}.temp.closest.tsv
		Rscript {params.codeDir}/TSS_close.R --variants {params.outDir}/{wildcards.GTExTissue}.temp.closest.tsv --TSSint {params.outDir}/{wildcards.GTExTissue}.temp.100kbInterval.tsv | gzip > {output.filteredVariantsProximal}
		
		rm {params.outDir}/{wildcards.GTExTissue}.temp.100kbInterval.tsv
		rm {params.outDir}/{wildcards.GTExTissue}.temp.closest.tsv
		""")
		
# filter enhancer predictions to invidual biosamples needed for predictions
rule filter_enhancer_predictions:
	input:
		predictions = lambda wildcards: config["predFile"][wildcards.method]
	params:
		outDir = config["outDir"]
	output:
		filteredPredictions = os.path.join(config["outDir"], "{method}", "{Biosample}.filteredPredictions.tsv.gz")
	run:
		shell(
		"""
		set +o pipefail;
		
		# filter predictions
		zcat {input.predictions} | csvtk cut -t -f chr,start,end,CellType,TargetGene | sed 1d | awk -v awkvar={wildcards.Biosample} '$4==awkvar' | cut -f 1,2,3,5 | sort -k1,1 -k2,2n | uniq | gzip > {output.filteredPredictions}
		""")

# generate list of variants with proximal genes and where the variant is located with respect to the enhancer predictions for each pairing of GTEx tissue/biosample 
rule make_prediction_table:
	input:
		filteredVariantsProximal = os.path.join(config["outDir"], "variantFilesForPrediction", "{GTExTissue}.filteredVariants.proximalGenes.tsv.gz"),
		filteredPredictions = lambda wildcards: os.path.join(config["outDir"], "{method}", "{Biosample}.filteredPredictions.tsv.gz")
	params:
		outDir = config["outDir"],
		codeDir = config["codeDir"]
	output:
		predTable = os.path.join(config["outDir"], "{method}", "{GTExTissue}.{Biosample}.predictionTable.tsv")
	run:
		shell(
		"""
		set +o pipefail;
		
		# intersect predictions with relevant variants 
		zcat {input.filteredVariantsProximal} | cut -f 1-4 | sort -k 1,1 -k2,2n | uniq > {params.outDir}/variantFilesForPrediction/{wildcards.GTExTissue}.temp.uniqueVariants.tsv
		zcat {input.filteredPredictions} | bedtools intersect -wa -wb -a {params.outDir}/variantFilesForPrediction/{wildcards.GTExTissue}.temp.uniqueVariants.tsv -b stdin | cut -f 4,8 > {params.outDir}/{wildcards.method}/temp.PredictionTargetGenes.tsv
		
		# classify variants
		zcat {input.filteredVariantsProximal} > {params.outDir}/variantFilesForPrediction/{wildcards.GTExTissue}.temp.allVariants.tsv
		
		Rscript {params.codeDir}/classify_enhancer_predictions.R --variants {params.outDir}/variantFilesForPrediction/{wildcards.GTExTissue}.temp.allVariants.tsv --pred {params.outDir}/{wildcards.method}/temp.PredictionTargetGenes.tsv > {output.predTable}
		
		rm {params.outDir}/variantFilesForPrediction/{wildcards.GTExTissue}.temp.uniqueVariants.tsv
		rm {params.outDir}/variantFilesForPrediction/{wildcards.GTExTissue}.temp.allVariants.tsv
		rm {params.outDir}/{wildcards.method}/temp.PredictionTargetGenes.tsv
		""")

# generate plot of prediction rates across all methods and tissues

predTables = []
for x in config["methods"]:
	# list of prediction tables
	predTables.extend(expand(os.path.join(config["outDir"], x, "{GTExTissue}.{Biosample}.predictionTable.tsv"), zip, GTExTissue=config["predRates"]["mapGTExTissues"][x], Biosample=config["predRates"]["mapBiosamples"][x]))
	
rule plot_prediction_rates:
	input:
		allTables=predTables
	params:
		codeDir = config["codeDir"]
	output:
		predictionPlot = os.path.join(config["outDir"], "predictionRates.pdf")
	run:
		shell("""
		set +o pipefail;
		
		Rscript {params.codeDir}/plot_prediction_rates.R --tables "{input.allTables}" --out {output.predictionPlot}
		""")