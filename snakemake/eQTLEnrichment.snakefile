# coding: utf-8

import pandas as pd
from os.path import join

pred_config = config["predictionsTable"]
preds_config_file = pd.read_table(pred_config).set_index("entry", drop=False)

# Gathering all the outputs for all sets of predictions and variants
rule all:
	input:
		expand(os.path.join(config["outDir"], "{pred}/test.txt"), pred=config["predictions"]),
		expand(os.path.join(config["outDir"], "{pred}/enhancerPredictions.sorted.bed"), pred=config["predictions"]),
		expand(os.path.join(config["outDir"], "{pred}/GTExVariants-enhancerPredictions.PIP0.5.distalNoncoding.tsv.gz"), pred=config["predictions"]),
		expand(os.path.join(config["outDir"], "{pred}/commonVariants-enhancerPredictions.distalNoncoding.tsv.gz"), pred=config["predictions"]),
		os.path.join(config["outDir"], "variantLists/GTEx.filtered.PIP0.5.tsv.gz"),
		os.path.join(config["outDir"], "variantLists/Partition.distalNoncoding.bed"),
		os.path.join(config["outDir"], "variantLists/GTEx.filtered.PIP0.5.distalNoncoding.tsv.gz"),
		os.path.join(config["outDir"], "variantLists/distalNoncoding.bg.SNPs.bed.gz"),
		os.path.join(config["outDir"], "variantLists/variantsPerGTExTissue.tsv"),
		expand(os.path.join(config["outDir"], "{pred}/predictionCellTypes.tsv"), pred=config["predictions"]),
		expand(os.path.join(config["outDir"], "{pred}/enrichmentTable.tsv"), pred=config["predictions"]),
		expand(os.path.join(config["outDir"], "{pred}/countMatrix.tsv"), pred=config["predictions"]),
		expand(os.path.join(config["outDir"], "{pred}/basesPerEnhancerSet.tsv"), pred=config["predictions"]),
		expand(os.path.join(config["outDir"], "{pred}/commonVarPerBiosample.tsv"), pred=config["predictions"]),
#		expand(os.path.join(config["outDir"], "{pred}/enrichmentHeatmap.full.pdf"), pred=config["predictions"]),
#		expand(os.path.join(config["outDir"], "{pred}/enrichmentHeatmap.aggregated.pdf"), pred=config["predictions"]),
#		os.path.join(config["outDir"], "/cdf.pdf"),
#		os.path.join(config["outDir"], "/density.pdf")


rule testRule:
	input:
		predFile = lambda wildcard: preds_config_file.loc[wildcard.pred, "predFile"]
	output:
		testOut = os.path.join(config["outDir"], "{pred}/test.txt")
	run:
		shell(
				"""
				echo {input.predFile} > {output.testOut}	
				set +o pipefail;
				zcat {input.predFile} | head >> {output.testOut}
				""")

rule sortPredictions:
	input:
		predFile = lambda wildcard: preds_config_file.loc[wildcard.pred, "predFile"]
	params:
		chrSizes = config["chrSizes"]
	output:
		predictionsSorted = os.path.join(config["outDir"], "{pred}/enhancerPredictions.sorted.bed"),
		samples = os.path.join(config["outDir"], "{pred}/predictionCellTypes.tsv")
	run:
		shell(
			"""
			set +o pipefail;
			
			# get list of samples
			zcat {input.predFile} | csvtk cut -t -f chr,start,end,CellType | sed 1d | awk '{{a[$4]++}} END{{for(i in a) printf i" "}}' > {output.samples}
			
			# sort predictions file
			zcat {input.predFile} | csvtk cut -t -f chr,start,end,CellType | sed 1d | bedtools sort -i stdin -faidx {params.chrSizes} > {output.predictionsSorted}
			""")
			
rule intersectPredictionsVariants:
	input:
		predictionsSorted = lambda wildcard: os.path.join(config["outDir"], preds_config_file.loc[wildcard.pred, "entry"], "enhancerPredictions.sorted.bed"),
		filteredGTExVarDistalNoncoding = os.path.join(config["outDir"], "variantLists/GTEx.filtered.PIP0.5.distalNoncoding.tsv.gz")
	output:
		predictionsVariantsInt = os.path.join(config["outDir"], "{pred}/GTExVariants-enhancerPredictions.PIP0.5.distalNoncoding.tsv.gz")
	run:
		shell(
			"""
			# columns of output: 1-3 (variant loc), 4 (variant hgID), 5 (variant tissue), 6 (variant gene ens id), 7 (variant PIP), 8-10 (enhancer loc), 11 (enhancer cell type)
			zcat {input.filteredGTExVarDistalNoncoding} | bedtools intersect -wa -wb -a stdin -b {input.predictionsSorted} | gzip > {output.predictionsVariantsInt}
			""")

rule filterGTExVariants:
	input:
		GTExVariants = config["GTExvariants"]
	params: 
		ABCgenes = config["ABCgenes"],
		commonVar = config["bgVariants"],
		partition = config["partition"],
		codeDir = config["codeDir"]
	output: 
		filteredGTExVar = os.path.join(config["outDir"], "variantLists/GTEx.filtered.PIP0.5.tsv.gz"),
	run:
		shell(
			"""
			set +o pipefail;
			# filter all variants by SUSIE, credible set, PIP0.5; print set of columns: 1-3 (loc), 4 (hgID), 5 (tissue), 6 (gene ens id), 7 (PIP)
			zcat {input.GTExVariants} | awk '$16>=0.5 && $17 != -1  && $9 == "SUSIE"' | awk '{{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $10 "\t" $11 "\t" $16;}}' | sort | uniq | Rscript {params.codeDir}filter_to_ABC_genes.R --genes {params.ABCgenes} --col 6| gzip > {output.filteredGTExVar}
			""")
			
rule getVarPerGTExTissue:
	input:
		filteredGTExVarDistalNoncoding = os.path.join(config["outDir"], "variantLists/GTEx.filtered.PIP0.5.distalNoncoding.tsv.gz")
	output: 
		variantsPerGTExTissue = os.path.join(config["outDir"], "variantLists/variantsPerGTExTissue.tsv")
	run:
		shell(
			"""
			set +o pipefail;
			# get counts per GTEx tissue
			zcat {input.filteredGTExVarDistalNoncoding} | sort -k5 | awk '{{a[$5]++}} END {{for(i in a) print i"\t"a[i]}}' > {output.variantsPerGTExTissue}
			""")
			
rule filterToDistalNoncoding:
	input:
		commonVar = config["bgVariants"],
		filteredGTExVar = os.path.join(config["outDir"], "variantLists/GTEx.filtered.PIP0.5.tsv.gz")
	params: 
		partition = config["partition"],
		codeDir = config["codeDir"]
	output: 
		partitionDistalNoncoding = os.path.join(config["outDir"], "variantLists/Partition.distalNoncoding.bed"),
		filteredGTExVarDistalNoncoding = os.path.join(config["outDir"], "variantLists/GTEx.filtered.PIP0.5.distalNoncoding.tsv.gz"),
		commonVarDistalNoncoding = os.path.join(config["outDir"], "variantLists/distalNoncoding.bg.SNPs.bed.gz"),
	run:
		shell(
			"""
		awk '$4=="ABC" || $4=="AllPeaks" || $4=="Other" || $4=="OtherIntron"' {params.partition} > {output.partitionDistalNoncoding}

		zcat {input.filteredGTExVar} | bedtools intersect -wa -a stdin -b {output.partitionDistalNoncoding} | gzip > {output.filteredGTExVarDistalNoncoding}
		cat {input.commonVar} | bedtools intersect -wa -a stdin -b {output.partitionDistalNoncoding} | gzip > {output.commonVarDistalNoncoding}
			""")
		
rule computeCountMatrix:
	input: 
		variantsPredictionsInt = lambda wildcard: os.path.join(config["outDir"],preds_config_file.loc[wildcard.pred, "entry"], "GTExVariants-enhancerPredictions.PIP0.5.distalNoncoding.tsv.gz"),
		predictionsSorted = lambda wildcard: os.path.join(config["outDir"], preds_config_file.loc[wildcard.pred, "entry"], "enhancerPredictions.sorted.bed"),
		samples = lambda wildcard: os.path.join(config["outDir"],preds_config_file.loc[wildcard.pred, "entry"], "predictionCellTypes.tsv")
	params:
		# filteredGTExVar = os.path.join(config["outDir"], "variantLists/GTEx.filtered.PIP0.5.distalNoncoding.tsv.gz"),
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

rule computeCommonVarOverlap:
	input:
		predictionsSorted = lambda wildcard: os.path.join(config["outDir"], preds_config_file.loc[wildcard.pred, "entry"], "enhancerPredictions.sorted.bed"),
		commonVarDistalNoncoding = os.path.join(config["outDir"], "variantLists/distalNoncoding.bg.SNPs.bed.gz"),
		samples = lambda wildcard: os.path.join(config["outDir"],preds_config_file.loc[wildcard.pred, "entry"], "predictionCellTypes.tsv")
	params:
	output:
		commonVarPerBiosample = os.path.join(config["outDir"], "{pred}/commonVarPerBiosample.tsv"),
		commonVarPredictionsInt = os.path.join(config["outDir"], "{pred}/commonVariants-enhancerPredictions.distalNoncoding.tsv.gz")
	run:
		shell(
			"""			
			set +o pipefail;
			
			# intersect common var with enhancer predictions
			# columns: 1-3 (variant loc), 4 (variant rsid), 5-7 (enhancer loc), 8 (enhancer cell type)
			zcat {input.commonVarDistalNoncoding} | bedtools intersect -wa -wb -a stdin -b {input.predictionsSorted} | gzip > {output.commonVarPredictionsInt}
			
			samples=$(cat {input.samples})
			# get common variants: biosample/celltype follwed by count
			for CellType in $samples
			do
				printf $CellType'\\t' >> {output.commonVarPerBiosample}
				zcat {output.commonVarPredictionsInt} | awk -v awkvar=$CellType '$8==awkvar' | wc -l >> {output.commonVarPerBiosample}
				printf '\\n' >> {output.commonVarPerBiosample}
			done
			""")

rule computeEnhancerSetSizes:
	input:
		predictionsSorted = lambda wildcard: os.path.join(config["outDir"], preds_config_file.loc[wildcard.pred, "entry"], "enhancerPredictions.sorted.bed"),
		samples = lambda wildcard: os.path.join(config["outDir"],preds_config_file.loc[wildcard.pred, "entry"], "predictionCellTypes.tsv")
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
				cat {input.predictionsSorted} | awk -v awkvar=$CellType '$4==awkvar' | sort -k1,1 -k2,2n | bedtools merge -i stdin | awk 'BEGIN {{FS=OFS="\\t"}} {{print $3-$2}}' | awk '{{s+=$1}}END{{print s}}' >> {output.basesPerEnhancerSet}
				printf "\\n" >> {output.basesPerEnhancerSet}
			done
			""")

rule computeEnrichmentMatrix:
	input: 
		countMatrix = lambda wildcard: directory(os.path.join(config["outDir"], preds_config_file.loc[wildcard.pred, "entry"],"countMatrix.tsv"))
	params:
		variantsPerGTExTissue = os.path.join(config["outDir"], "variantLists/variantsPerGTExTissue.tsv"),
		commonVarDistalNoncoding = os.path.join(config["outDir"], "variantLists/distalNoncoding.bg.SNPs.bed.gz"),
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
			Rscript {params.codeDir}counts_to_enrichment.R --counts {input.countMatrix} --commonvar {params.commonVarPerBiosample} --varbytissue {params.variantsPerGTExTissue} --totalcommonvar $totalCommonVar > {output.enrichmentTable}
			""")
			
rule plotEnrichmentHeatmaps:
	input: 
		enrichmentTable = lambda wildcard: os.path.join(config["outDir"], preds_config_file.loc[wildcard.pred, "entry"],"enrichmentTable.tsv")
	params:
		sampleKey = lambda wildcard: preds_config_file.loc[wildcard.pred, "sampleKey"],
		sampleID = lambda wildcard: preds_config_file.loc[wildcard.pred, "sampleID"],
		sampleName = lambda wildcard: preds_config_file.loc[wildcard.pred, "sampleName"],
		sampleCategory = lambda wildcard: preds_config_file.loc[wildcard.pred, "sampleCategory"],
		basesPerEnhancerSet = os.path.join(config["outDir"], "{pred}/basesPerEnhancerSet.tsv"),
		codeDir = config["codeDir"],
		
	output: 
		heatmapFull = os.path.join(config["outDir"], "{pred}/enrichmentHeatmap.full.pdf"),
		heatmapAggregated = os.path.join(config["outDir"], "{pred}/enrichmentHeatmap.aggregated.pdf")

	run:
		shell(
			"""	
			set +o pipefail;
			
			# plot full heat map
			Rscript {params.codeDir}plot_enrichment_heatmap.R --table {input.enrichmentTable} --outfile {output.heatmapFull} --samplekey {params.sampleKey} --cellid {params.sampleID} --category {params.sampleName} --enhancersizes {params.basesPerEnhancerSet}
			
			# plot aggregated heat map
			Rscript {params.codeDir}plot_enrichment_heatmap.R --table {input.enrichmentTable} --outfile {output.heatmapAggregated} --samplekey {params.sampleKey} --cellid {params.sampleID} --category {params.sampleCategory} --enhancersizes {params.basesPerEnhancerSet}
			
			""")
			
rule plotComparisons:
	input: 
		enrichmentTables = expand(os.path.join(config["outDir"],"{pred}/enrichmentTable.tsv"), pred=config["predictions"]),
	params:
		names = config["predictions"],
		codeDir = config["codeDir"],
		outDir = config["outDir"]
	output: 
		cdf = os.path.join(config["outDir"], "/cdf.pdf"),
		density = os.path.join(config["outDir"], "/density.pdf")
	run:
		shell(
			"""	
				set +o pipefail;
				Rscript {params.codeDir}plot_comparisons.R --names {params.names} --tables {input.enrichmentTables} --outdir {params.outDir}
			
			""")
			
