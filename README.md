# eQTL benchmarking pipeline
We developed this "eQTL benchmarking pipeline” to assess whether predictive models of enhancer-gene regulatory interactions could (i) identify enhancers enriched for fine-mapped eQTL variants, and (ii) link those enhancers to their corresponding eQTL target genes. We analyzed variants from the GTEx V8 resource26 fine-mapped using the Sum of Single Effects (SuSIE) model.

For this analysis, we focus on distal noncoding eQTLs by by filtering out (i) coding sequences, 5′ and 3′ untranslated regions of protein-coding genes, and splice sites (within 10 bp of a intron–exon junction of a protein-coding gene) of protein-coding genes, and (ii) promoters (±250 bp from the gene TSS) of protein-coding genes. We define **enrichment** as (fraction of eQTL variants with PIP >= 50% that overlap predicted enhancers) / (fraction of all 1000G SNPs that overlap predicted enhancers). We define **recall (linking)** as as the fraction of variants overlapping a predicted enhancer linked to at least one of its eGenes, and **recall (total)** as the fraction of variants overlapping predicted enhancers

## Outputs
Running this pipeline will produce:
1. Comparison of the enrichment of eQTL variants versus common variants, recall (not linking), and fraction of variants overlapping a predicted enhancer linked to the variant's eGene given that the variant overlaps an enhancer. Here, predicted enhancers are defined as enhancer-gene links with a score greater than the provided threshold value. For each perdictor, each point represents variants from one GTEx tissue and predicted enhancers from the specified matching biosample. These metrics will be calculated and displayed at the specified distance stratifications of eGene-eVariant separation.
2. Enrichment-recall curves for each GTEx tissue with specified biosample pairings and for all tissue/biosample matches, where enrichment and recall (linking) are calculated at the specified number of threshold values yielding even steps across the range of recall for that method
3. Comparison of enrichment at specified recall values for each GTEx tissue with specified pairings across predictors
4. Heatmaps showing enrichment at every tissue/biosample pairing for each method

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
	- **eQTLVariants:** gzipped, tab-separated file of eQTL variants with the following columns: `chr`, `start`, `end`, `varID_hg38` (variant ID), `gene_hgnc` (eGene), `tissue`, `pip`. Here, we use fine-mapped variants from GTEx that have been filtered for variants in a credible set and expressed with a TPM>1 in their respective tissue
		- The complete file of fine-mapped GTEx variants is available on Synapse [here](https://www.synapse.org/#!Synapse:syn52264297), and the expression data from GTEx is available [here](https://www.synapse.org/#!Synapse:syn52264240). The code to process these files into the format for this pipeline is under `workflow/scripts/generate_inputs`.
	- Genomic annotations **chrSizes**, **partition**, **TSS**, **GTExGeneUniverse** (all provided here in the `resources` directory)
	- **bgVariants:** a .bed file of common SNPs  with columns (no header) `chr`, `start`, `end`, `rsid`. The list of background variants we use can be downloaded on Synapse [here](https://www.synapse.org/#!Synapse:syn52264319) and was collated from [here](https://alkesgroup.broadinstitute.org/LDSCORE/baseline_v1.1_hg38_annots/).
	- The following parameters:
		- **distances:** a list of values, in Kb, for eGene-eVariant distance-stratified benchmarking
		- **recalls:** a list of recall values at which to compare enrichments
		- **thresholdPIP:** PIP threshold to filter variants
		- **nThresholdSteps:** the number of thresholds to include when generating enrichment-recall curves. We recommend using about 20 steps.
  
2. **Methods config table** - an example file is included at: `config/methods_config_example.tsv`. The following columns (with a header) are required, where each row represents one predictive method:
   - **method:** predictive method name, to match listed methods in the main config file outlined above
   - **sampleKey:** path to .tsv file with columns and header `biosample`, `predictionFile`, and `GTExTissue`. The `GTExTissue` column should be left EMPTY for biosamples without a match. Multiple GTEx tissues can be matched with a single biosample in a comma-separated list. The prediction files must be separated by biosample and must have a header with minimally the columns `chr`, `start`, `end`, `TargetGene`, `score_col` (see below)
   - **geneUniverse:** path to .bed6 file specifying which genes are considered by this method. Variants will be filtered to this set of genes for the analysis.
   - **pred_name_long:** a string with the full method name to be used in plots
   - **threshold:** score threshold value to be used to generate boxplots stratified by distance and enrichment heat maps (outputs (1) and (4)). We used the threshold values corresponding to 70% recall from our CRISPR benchmarking pipeline for each model.
   - **score_col:** name of the column in prediction files with the prediction score
   - **color:** hex code specifying what color to use in plots for this method. May be left blank and will be automatically filled.
   - **inverse_predictor:** `TRUE` if high scores correspond to lower prediction confidence (eg. distance to TSS), otherwise `FALSE`
   - **boolean:** `TRUE` if this is a binary 0 or 1 predictor, otherwise `FALSE`


