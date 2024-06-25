# eQTL bechmarking pipeline
We developed this "eQTL benchmarking pipeline” to assess whether predictive models of enhancer-gene regulatory interactions could (i) identify enhancers enriched for fine-mapped eQTL variants, and (ii) link those enhancers to their corresponding eQTL target genes. We have applied this pipeline using eQTLs fine-mapped with the Sum of Single Effects (SuSiE) model from the [GTEx V8 resource](https://gtexportal.org/home/downloads/adult-gtex/overview)  and from the [eQTL Catalogue release 7](https://www.ebi.ac.uk/eqtl/Data_access/).

For this analysis, we focus on distal noncoding eQTLs by by filtering out (i) coding sequences, 5′ and 3′ untranslated regions of protein-coding genes, and splice sites (within 10 bp of a intron–exon junction of a protein-coding gene) of protein-coding genes, and (ii) promoters (±250 bp from the gene TSS) of protein-coding genes. We define **enrichment** as (fraction of eQTL variants with PIP >= 50% that overlap predicted enhancers) / (fraction of all 1000G SNPs that overlap predicted enhancers). We define **recall (linking)** as as the fraction of eVariant-eGene pairs overlapping a predicted enhancer linked to the correct eGenes, and **recall (total)** as the fraction of eVariant-eGene pairs with the eVariant overlapping predicted a enhancer.

## Outputs
Running this pipeline will produce the following plots and corresponding tables, many of which are also included in a summary HTML report (`results/benchmarkingReport.html`). 
1. To compare comprehensively compare performance across predictors, we plot the enrichment and recall for user-specified eQTL biosample, prediction biosample matches, where predicted enhancers are defined as enhancer-gene links with a score greater than the provided threshold value.
	- `results/plots/thresholdedPerformanceComparison.pdf`: Boxplots for for enrichment, recall (total), and fraction of variants overlapping a predicted enhancer linked to the variant's eGene given that the variant overlaps an enhancer for eQTLs stratified by eGene-eVariant distance. Each datapoint represents variants from one eQTL biosample and predicted enhancers from the specified prediction biosample. Metrics for for eQTL biosamples with fewer than 20 variants are removed, and outliers are plotted as points.
	- `results/plots/allMatchedThresholdedPerformance.scatter.pdf`: Scatter plot of enrichment and recall (linking) for matched eQTL and prediction biosamples for variants with any eVariant-eGene distance. Points for eQTL biosamples with fewer than 20 variants are removed as above is performed.
2. To evaluate cell-type specificity of predictors, we plot enrichment and recall (linking) heatmaps across all eQTL biosamples, prediction biosample pairs, where predicted enhancers are defined as above based on the provided threshold value. Each axis is clustered based on Pearson correlation The number of base pairs predicted to be in enhancers for each biosample and the number of variants (enrichment) or eVariant-eGene pairs (recall) per eQTL biosample are plotted as well.
	- `results/plots/enrichmentHeatmaps/{predictor}.enrichmentHeatmap.withMetrics.pdf`: Color scale is dynamic based on the range of enrichments for that predictor. Stars are plotted at biosample intersections with an enrichment significantly greater than 1 (hypergeometric test, Bonferroni-adjusted p-value).
	- `results/plots/recallHeatmaps/{predictor}.recallHeatmap.withMetrics.pdf`:  Color scale is constant across all predictors from 0 to 0.25. eQTL biosamples with constant recall for all prediction biosamples are removed. 
3. To evaluate predictors agnostic to a predefined threshold and account for the tradeoff between high enrichment and low recall, we calculate enrichment and recall (linking) curves across the full range of threshold values for specified eQTL biosample, prediction biosample matches. 
	- `results/plots.enrichmentRecall/enrichmentRecall.GTExTissue{tissue}.pdf`: One plot is created for each eQTL biosample with at least one specified prediction biosample match. One curve is plotted per prediction biosample matching the eQTL biosample, with straight lines connecting datapoints evaluated at thresholds yielding even steps across the recall range for all specified biosample matches. If multiple biosamples for the same predictor match one eQTL biosample, the curves will vary slightly in shade and the key will specify the specific match. Binary predictors are plotted as points corresponding to the positive predicted enhancers. The enrichment axis is cut off at 50 to eliminate spiking at very high thresholds.
	- `results/plots.enrichmentRecall/enrichmentRecall.GTExTissueAllMatches.pdf`: Enrichment-recall curve as above where values are aggregated across all specified biosample matches for each predictor. Enrichment and recall are both aggregated by summing variant overlap counts across each biosample match before computing the metric.
	- `results/plots/enrichmentAtRecall/enrichmentAtRecall{recall}.GTExTissue{tissue}.pdf`: For user-specified recall (linking) values, a plot is generated for each eQTL biosample with at least one predictor biosample designated match. Each predictor biosample match is only included in the plot if they acheive a recall (linking) within 0.02 of the desired value; the exact recall corresponding to the threshold at which enrichment is plotted is shown in the legend. We calculate Bonferroni-adjusted p-values between all pairwise comparisons of included predictors using the approach of [comparing relative risks](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1125071/). Pairs with significant differences are not indicated on the plot because algorithimically calculating the y-position of the required bars is confusing. If fewer than two predictor biosample matches meet the designated recall, empty placeholder files are created. 

## Running the pipeline
1. Clone this repository
```
git@github.com:EngreitzLab/eQTLEnrichment.git
```
2. Edit configuration files (see below)
3. Activate a conda environment with [mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html) and [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) installed.
4. Run pipeline
```
snakemake -j1 --configfile config/config_example.yml --use-conda
```

## Config files
The pipeline requires two config files. The required inputs are outlined below.
1. **Main config file** - an example file is included at: `config/config_example.yml`. The following fields are required:
	- **outDir:** output directory name
	- **methodsTable:** configuration file with information about predictors (see below for specifications)
	- **methods:** list of methods from `methodsTable` to be benchmarked
	- **eQTLVariants:** gzipped, tab-separated file of eQTL variants with the following columns: `chr`, `start`, `end`, `varID_hg38` (variant ID), `gene_hgnc` (eGene), `tissue`, `pip`. Here, we use fine-mapped variants from GTEx or the eQTL Catalogue. For the GTEx dataset, we implement additional filtering for variants in a credible set and expressed with a TPM>1 in their respective tissue 
		- The unprocessed file of SuSiE-fine-mapped GTEx variants is available on Synapse [here](https://www.synapse.org/#!Synapse:syn52264297), and the expression data from GTEx is available [here](https://www.synapse.org/#!Synapse:syn52264240). The code to process these files into the format for this pipeline is under `workflow/scripts/generate_inputs`.
		- The code to generate an input file from the eQTL Catalogue is available upon request.
	- Genomic annotations **chrSizes**, **partition**, **TSS**, **GTExGeneUniverse** (all provided here in the `resources` directory)
	- **bgVariants:** a .bed file of 1000G SNPs  with columns (no header) `chr`, `start`, `end`, `rsid`. The list of background variants we use can be downloaded on Synapse [here](https://www.synapse.org/#!Synapse:syn52264319) and was collated from https://alkesgroup.broadinstitute.org/LDSCORE/baseline_v1.1_hg38_annots/
	- The following parameters:
		- **distances:** a list of values, in Kb, for eGene-eVariant distance-stratified benchmarking
		- **recalls:** a list of recall values at which to compare enrichments. May be left blank.
		- **thresholdPIP:** PIP threshold to filter variants
		- **thresholdPval:"** p-value threshold used in statistical tests
		- **nThresholdSteps:** number of thresholds to include when generating enrichment-recall curves. We recommend using 20-50 steps.
  
2. **Methods config table** - an example file is included at: `config/methods_config_example.tsv`. The following columns (with a header) are required, where each row represents one predictive method:
   - **method:** predictive method name, to match listed methods in the main config file outlined above
   - **sampleKey:** path to .tsv file with columns and header `biosample`, `predictionFile`, and `GTExTissue`. The `GTExTissue` column should be left EMPTY for biosamples without a match. Multiple GTEx tissues can be matched with a single biosample in a comma-separated list. The prediction files must be separated by biosample and must have a header with minimally the columns `chr`, `start`, `end`, `TargetGene`, and a numeric score column (see below)
   - **geneUniverse:** path to .bed6 file specifying which genes are considered by this method. Variants will be filtered to this set of genes for the analysis.
   - **pred_name_long:** a string with the full method name to be used in plots
   - **threshold:** score threshold value to be used to generate boxplots stratified by distance and enrichment heat maps (outputs (1) and (4)). We used the threshold values corresponding to 70% recall from our CRISPR benchmarking pipeline for each model.
   - **score_col:** name of the column in prediction files with the prediction score
   - **color:** hex code specifying what color to use in plots for this method. May be left blank and will be automatically filled.
   - **inverse_predictor:** `TRUE` if high scores correspond to lower prediction confidence (eg. distance to TSS), otherwise `FALSE`
   - **boolean:** `TRUE` if this is a binary 0 or 1 predictor, otherwise `FALSE`

## Works-in-progress
- Edit terminology in configuration, code, file names to not be specific to GTEx (replace "GTEx tissue" with "eQTL biosample")
- Implement benchmarking of groups of prediction biosamples against a single eQTL biosamples to better map to hetergenous tissues
- Integrate computations from enrichment and recall into a single function
- Implement baseline predictors computed internally (e.g. distance to TSS, random expressed gene with 500kb)
- Review/improve the "enrichment at recall" analysis
