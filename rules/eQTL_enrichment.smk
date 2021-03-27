# sort enhancer predictions by chromosome & start location, get list of samples
rule sortPredictions:
	input:
		predFile = lambda wildcards: config["predFile"][wildcards.method]
	params:
		chrSizes = config["chrSizes"]
	output:
		predictionsSorted = os.path.join(config["outDir"], "{method}", "enhancerPredictions.sorted.bed.gz"),
		samples = os.path.join(config["outDir"], "{method}", "biosampleList.tsv")
	run:
		shell(
			"""
			set +o pipefail;
			
			# get list of samples
			zcat {input.predFile} | csvtk cut -t -f chr,start,end,CellType | sed 1d | awk '{{a[$4]++}} END{{for(i in a) printf i"\\t"}}' > {output.samples}
			
			# sort predictions file
			zcat {input.predFile} | csvtk cut -t -f chr,start,end,CellType | sed 1d | sort -k1,1 -k2,2n | gzip > {output.predictionsSorted}
			""")

# intersect predictions and GTEx variants			
rule intersectPredictionsVariants:
	input:
		predictionsSorted = os.path.join(config["outDir"], "{pred}", "enhancerPredictions.sorted.bed.gz"),
		filteredGTExVarDistalNoncoding = os.path.join(config["outDir"], "variantFilesForEnrichment/GTEx.filtered.PIP0.5.distalNoncoding.tsv.gz")
	params: 
		outDir = config["outDir"]
	output:
		predictionsVariantsInt = os.path.join(config["outDir"], "{pred}/GTExVariants-enhancerPredictions.PIP0.5.distalNoncoding.tsv.gz")
	run:
		shell(
			"""
			# columns of output: 1-3 (variant loc), 4 (variant hgID), 5 (variant tissue), 6 (variant gene ens id), 7 (variant PIP), 8-10 (enhancer loc), 11 (enhancer cell type)
			zcat {input.filteredGTExVarDistalNoncoding} > {params.outDir}temp.tsv
			zcat {input.predictionsSorted} | bedtools intersect -wa -wb -sorted -a {params.outDir}temp.tsv -b stdin | gzip > {output.predictionsVariantsInt}
			rm {params.outDir}temp.tsv
			
			""")

# filter GTEx variants by PIP, credible set, protein-coding genes
rule filterGTExVariants:
	input:
		GTExVariants = config["GTExVariants"]
	params: 
		ABCgenes = config["genes"],
		commonVar = config["bgVariants"],
		partition = config["partition"],
		codeDir = config["codeDir"]
	output: 
		filteredGTExVar = os.path.join(config["outDir"], "variantFilesForEnrichment/GTEx.filtered.PIP0.5.tsv.gz"),
	run:
		shell(
			"""
			set +o pipefail;
			# filter all variants by SUSIE, credible set, PIP0.5; print set of columns: 1-3 (loc), 4 (hgID), 5 (tissue), 6 (gene ens id), 7 (PIP)
			zcat {input.GTExVariants} | awk '$16>=0.5 && $17 != -1  && $9 == "SUSIE"' | awk '{{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $10 "\t" $11 "\t" $16;}}' | sort -k1,1 -k2,2n | uniq | Rscript {params.codeDir}/filter_to_ABC_genes.R --genes {params.ABCgenes} --col 6| gzip > {output.filteredGTExVar}
			""")

# get number of variants per GTEx tissue
rule getVarPerGTExTissue:
	input:
		filteredGTExVarDistalNoncoding = os.path.join(config["outDir"], "variantFilesForEnrichment/GTEx.filtered.PIP0.5.distalNoncoding.tsv.gz")
	output: 
		variantsPerGTExTissue = os.path.join(config["outDir"], "variantFilesForEnrichment/variantsPerGTExTissue.tsv")
	run:
		shell(
			"""
			set +o pipefail;
			# get counts per GTEx tissue
			zcat {input.filteredGTExVarDistalNoncoding} | sort -k5 | awk '{{a[$5]++}} END {{for(i in a) print i"\t"a[i]}}' > {output.variantsPerGTExTissue}
			""")

# filter GTEx variants and common variants to distal noncoding
rule filterToDistalNoncoding:
	input:
		commonVar = config["bgVariants"],
		filteredGTExVar = os.path.join(config["outDir"], "variantFilesForEnrichment/GTEx.filtered.PIP0.5.tsv.gz")
	params: 
		partition = config["partition"],
		codeDir = config["codeDir"]
	output: 
		partitionDistalNoncoding = os.path.join(config["outDir"], "variantFilesForEnrichment/Partition.distalNoncoding.bed"),
		filteredGTExVarDistalNoncoding = os.path.join(config["outDir"], "variantFilesForEnrichment/GTEx.filtered.PIP0.5.distalNoncoding.tsv.gz"),
		commonVarDistalNoncoding = os.path.join(config["outDir"], "variantFilesForEnrichment/distalNoncoding.bg.SNPs.bed.gz"),
	run:
		shell(
			"""
		awk '$4=="ABC" || $4=="AllPeaks" || $4=="Other" || $4=="OtherIntron"' {params.partition} | sort -k1,1 -k2,2n > {output.partitionDistalNoncoding}

		zcat {input.filteredGTExVar} | bedtools intersect -wa -sorted -a stdin -b {output.partitionDistalNoncoding} | gzip > {output.filteredGTExVarDistalNoncoding}
		cat {input.commonVar} | bedtools intersect -wa -sorted -a stdin -b {output.partitionDistalNoncoding} | gzip > {output.commonVarDistalNoncoding}
			""")
		
# compute matrix with number of variants from each GTEx tissue overlapping enhancers from each biosample
rule computeCountMatrix:
	input: 
		#variantsPredictionsInt = lambda wildcard: os.path.join(config["outDir"],preds_config_file.loc[wildcard.pred, "entry"], "GTExVariants-enhancerPredictions.PIP0.5.distalNoncoding.tsv.gz"),
		#predictionsSorted = lambda wildcard: os.path.join(config["outDir"], preds_config_file.loc[wildcard.pred, "entry"], "enhancerPredictions.sorted.bed.gz"),
		#samples = lambda wildcard: os.path.join(config["outDir"],preds_config_file.loc[wildcard.pred, "entry"], "predictionCellTypes.tsv")
		variantsPredictionsInt = os.path.join(config["outDir"],"{pred}", "GTExVariants-enhancerPredictions.PIP0.5.distalNoncoding.tsv.gz"),
		predictionsSorted = os.path.join(config["outDir"], "{pred}", "enhancerPredictions.sorted.bed.gz"),
		samples = os.path.join(config["outDir"], "{pred}", "biosampleList.tsv")
	params:
		# filteredGTExVar = os.path.join(config["outDir"], "variantFilesForEnrichment/GTEx.filtered.PIP0.5.distalNoncoding.tsv.gz"),
		outDir = config["outDir"]
	output: 
		countMatrix = os.path.join(config["outDir"], "{pred}/countMatrix.tsv"),
	run:
		shell(
			"""			
			set +o pipefail;
			
			samples=$(cat {input.samples})
			for CellType in $samples
			do
				printf $CellType"\\t" >> {output.countMatrix}
				# columns of intersection file: 1-3 (variant loc), 4 (variant hgID), 5 (variant tissue), 6 (variant gene ens id), 7 (variant PIP), 8-10 (enhancer loc), 11 (enhancer cell type)
				# filter for cell type, select just variant hgID and tissue, deduplicate (through uniq)
				# print counts per GTEx tissue, sort (alphabetically by first tissue), print only count, then print all numbers in a row in matrix
				zcat {input.variantsPredictionsInt} | awk -v awkvar=$CellType '$11==awkvar' | cut -f 4,5 | sort -k1 | uniq | awk '{{a[$2]++}} END {{for(i in a) print i" "a[i]}}' | sort -k1 | cut -d' ' -f2 | awk 'BEGIN {{ ORS = "\\t" }} {{ print }}' >> {output.countMatrix}
				printf "\\n" >> {output.countMatrix}

			done
			""")
			
# computer number of common variants overlapping enhancers in each biosample
rule computeCommonVarOverlap:
	input:
		# predictionsSorted = lambda wildcard: os.path.join(config["outDir"], preds_config_file.loc[wildcard.pred, "entry"], "enhancerPredictions.sorted.bed.gz"),
		predictionsSorted = os.path.join(config["outDir"], "{pred}", "enhancerPredictions.sorted.bed.gz"),
		commonVarDistalNoncoding = os.path.join(config["outDir"], "variantFilesForEnrichment/distalNoncoding.bg.SNPs.bed.gz"),
		# samples = lambda wildcard: os.path.join(config["outDir"],preds_config_file.loc[wildcard.pred, "entry"], "predictionCellTypes.tsv")
		samples = os.path.join(config["outDir"], "{pred}", "biosampleList.tsv")
		
	params:
		outDir = config["outDir"]
	output:
		commonVarPerBiosample = os.path.join(config["outDir"], "{pred}/commonVarPerBiosample.tsv"),
		commonVarPredictionsInt = os.path.join(config["outDir"], "{pred}/commonVariants-enhancerPredictions.distalNoncoding.tsv.gz")
	run:
		shell(
			"""			
			set +o pipefail;
			
			# intersect common var with enhancer predictions
			# columns: 1-3 (variant loc), 4 (variant rsid), 5-7 (enhancer loc), 8 (enhancer cell type)
			zcat {input.commonVarDistalNoncoding} > {params.outDir}temp.tsv
			zcat {input.predictionsSorted} | bedtools intersect -wa -wb -sorted -a {params.outDir}temp.tsv -b stdin | gzip > {output.commonVarPredictionsInt}
			rm {params.outDir}temp.tsv
			
			samples=$(cat {input.samples})
			# get common variants: biosample/celltype followed by count
			for CellType in $samples
			do
				printf $CellType'\\t' >> {output.commonVarPerBiosample}
				zcat {output.commonVarPredictionsInt} | awk -v awkvar=$CellType '$8==awkvar' | wc -l >> {output.commonVarPerBiosample}
				# printf '\\n' >> {output.commonVarPerBiosample}
			done
			""")

# compute number of base pairs in each enhancer set
rule computeEnhancerSetSizes:
	input:
		#predictionsSorted = lambda wildcard: os.path.join(config["outDir"], preds_config_file.loc[wildcard.pred, "entry"], "enhancerPredictions.sorted.bed.gz"),
		#samples = lambda wildcard: os.path.join(config["outDir"],preds_config_file.loc[wildcard.pred, "entry"], "predictionCellTypes.tsv")
		predictionsSorted = os.path.join(config["outDir"], "{pred}", "enhancerPredictions.sorted.bed.gz"),
		samples = os.path.join(config["outDir"],"{pred}", "biosampleList.tsv")
	params:
	output:
		basesPerEnhancerSet = os.path.join(config["outDir"], "{pred}/basesPerEnhancerSet.tsv"),
	run:
		shell(
			"""			
			set +o pipefail;
			
			samples=$(cat {input.samples})
			
			for CellType in $samples
			do
				# make list of sizes of enhancer sets for each biosample
				printf $CellType"\\t" >> {output.basesPerEnhancerSet}
				zcat {input.predictionsSorted} | awk -v awkvar=$CellType '$4==awkvar' | sort -k1,1 -k2,2n | bedtools merge -i stdin | awk 'BEGIN {{FS=OFS="\\t"}} {{print $3-$2}}' | awk '{{s+=$1}}END{{print s}}' >> {output.basesPerEnhancerSet}
				# printf "\\n" >> {output.basesPerEnhancerSet}
			done
			""")
			
# generate matrix with enrichment values for each GTEx tissue/biosample intersection
rule computeEnrichmentMatrix:
	input: 
		# countMatrix = lambda wildcard: directory(os.path.join(config["outDir"], preds_config_file.loc[wildcard.pred, "entry"],"countMatrix.tsv"))
		countMatrix = os.path.join(config["outDir"], "{pred}","countMatrix.tsv")
	params:
		variantsPerGTExTissue = os.path.join(config["outDir"], "variantFilesForEnrichment/variantsPerGTExTissue.tsv"),
		commonVarDistalNoncoding = os.path.join(config["outDir"], "variantFilesForEnrichment/distalNoncoding.bg.SNPs.bed.gz"),
		commonVarPerBiosample = os.path.join(config["outDir"], "{pred}/commonVarPerBiosample.tsv"),
		basesPerEnhancerSet = os.path.join(config["outDir"], "{pred}/basesPerEnhancerSet.tsv"),
		codeDir = config["codeDir"],
		outDir = config["outDir"]
	output: 
		enrichmentTable = os.path.join(config["outDir"], "{pred}/enrichmentTable.tsv")
	run:
		shell(
			"""			
			set +o pipefail;

			# get total common variants
			totalCommonVar=$(zcat {params.commonVarDistalNoncoding} | wc -l)

			# pass to R to calc and format enrichment matrix
			Rscript {params.codeDir}/counts_to_enrichment.R --counts {input.countMatrix} --commonvar {params.commonVarPerBiosample} --varbytissue {params.variantsPerGTExTissue} --totalcommonvar $totalCommonVar > {output.enrichmentTable}
			""")

# plot aggregate and full heat maps of enrichment values
rule plotEnrichmentHeatmaps:
	input: 
		enrichmentTable = os.path.join(config["outDir"], "{pred}","enrichmentTable.tsv")
		
	params:
		basesPerEnhancerSet = os.path.join(config["outDir"], "{pred}/basesPerEnhancerSet.tsv"),
		codeDir = config["codeDir"],
		sampleKey = lambda wildcards: config["enrichment_heatmaps"]["sampleKey"][wildcards.pred],
		sampleID = lambda wildcards: config["enrichment_heatmaps"]["sampleID"][wildcards.pred],
		sampleName = lambda wildcards: config["enrichment_heatmaps"]["sampleName"][wildcards.pred],
		sampleCategory = lambda wildcards: config["enrichment_heatmaps"]["sampleCategory"][wildcards.pred]
		
	output: 
		heatmapFull = os.path.join(config["outDir"], "{pred}/enrichmentHeatmap.full.pdf"),
		heatmapAggregated = os.path.join(config["outDir"], "{pred}/enrichmentHeatmap.aggregated.pdf")

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
rule plotComparisons:
	input: 
		enrichmentTables = expand(os.path.join(config["outDir"],"{pred}/enrichmentTable.tsv"), pred=config["methods"]),
	params:
		names = config["methods"],
		codeDir = config["codeDir"],
		outDir = config["outDir"]
	output: 
		cdf = os.path.join(config["outDir"], "cdf.pdf"),
		density = os.path.join(config["outDir"], "density.pdf")
	run:
		shell(
			"""	
				set +o pipefail;
				Rscript {params.codeDir}/plot_comparisons.R --names "{params.names}" --tables "{input.enrichmentTables}" --outdir {params.outDir}
			
			""")
			
