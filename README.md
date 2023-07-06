# eQTL bechmarking pipeline
We developed this "eQTL benchmarking pipeline” to assess whether predictive models of enhancer-gene regulatory interactions could (i) identify enhancers enriched for fine-mapped eQTL variants, and (ii) link those enhancers to their corresponding eQTL target genes. We analyzed variants from the GTEx V8 resource26 fine-mapped using the Sum of Single Effects (SuSIE) model.

For this analysis, we focus on distal noncoding eQTLs by by filtering out (i) coding sequences, 5′ and 3′ untranslated regions of protein-coding genes, and splice sites (within 10 bp of a intron–exon junction of a protein-coding gene) of protein-coding genes, and (ii) promoters (±250 bp from the gene TSS) of protein-coding genes. We define **enrichment** as (fraction of eQTL variants with PIP >= 50% that overlap predicted enhancers) / (fraction of all 1000G SNPs that overlap predicted enhancers). We define **recall (linking)** as as the fraction of variants overlapping a predicted enhancer linked to at least one of its eGenes, and **recall (not linking)** as the fraction of variants overlapping predicted enhancers

## Outputs
Running this pipeline will produce:
1. Boxplots showing the enrichment of eQTL variants versus common variants, recall (not linking), and fraction of variants overlapping a predicted enhancer linked to the variant's eGene given that the variant overlaps an enhancer. Predicted enhancers are defined as predictions with a score greater than the specified threshold value. For each method, each point represents variants from one GTEx tissue and predicted enhancers from the specified matching biosample. These metrics will be calculated and displayed at the specified distance stratifications of eGene-eVariant separation.
2. Enrichment-recall curves for all specified GTEx tissue-biosample pairings, where enrichment and recall (linking) are calculated at the specified number of threshold values yielding even steps across the range of recall for that method.
3. Heatmaps showing enrichment at every tissue/biosample pairing for each method

**Notes:** Outputs (1) and (2) are saved as individual pdfs in  `/plots` of your output directory and can also be viewed in the final eQTL benchmarking report. (1) is only produced if there is more than one GTEx tissue paired with a biosample across the methods. If multiple biosamples are matched with the same tissue, only one (2) will be produced (from the first biosample listed in the `sampleKey` matching that tissue). Output (3) can be viewed in the benchmarking report.

## Running the pipeline
To run the pipeline, call each of the rules `first` through `fifth` in consecutive order as shown in the example `workflow/run-GM12878_predictors.sh` file. Each rule is dependent on successful completion of the previous ones to successfully run. Note that at this time, rule `fifth`, which generates the final benchmarking report, will not run if there is only one tissue/biosample match.

## Config files
The pipeline requires two config files. The required inputs are outlined below.
1. **Main config file** - an example file is included at: `config/config-hg38-GM12878_predictors.yml`. The following fields are required:
    - **outDir:** output directory name
    - **codeDir:** directory with scripts for running the pipeline (should point to `workflow/scripts/`)
    - **envDir:** directory with conda environments (should point to `workflow/envs`)
    - **methodsTable:** path to methods config file (see below for specifications)
    - **GTExVariants:** path to file of fine-mapped GTEx variants, with following columns and no header. If using a file with another format, edit the column numbers and filtering strategy specified in the rule `filter_all_variants` in `workflow/rules/preprocessing.smk`. The file of fine-mapped GTEx variants used in the ENCODE analysis can be downloaded here: https://mitra.stanford.edu/kundaje/kmualim/out/Benchmarking/GTEx/GTEx_30tissues_hg38.tsv.gz or on Synapse here: https://doi.org/10.7303/syn18134906. 
      - *chromosome:*      chromosome (autosomes only)
      - *start:*   start position of variant (0-indexed)
      - *end:*     end position of variant (0-indexed)
      - *variant:* unique variant identifier in hg19 (chr:pos:ref:alt)
      - *variant_hg38:*    unique variant identifier in hg38 (chr:pos:ref:alt)
      - *allele1:* reference allele 
      - *allele2:* alternate allele
      - *cohort:*  GWAS cohort
      - *method:*  fine-mapping method used
      - *tissue:*  name of tissue where gene expression was measured
      - *gene:*    name of gene of interest for cis-eQTLs (Ensembl ID)
      - *maf:*     allele frequency of the minor allele in cohort
      - *beta_marginal:*   marginal association effect size from linear mixed model
      - *se_marginal: *    standard error on marginal association effect size from linear mixed model
      - *z:*       test statistic for marginal association
      - *pip:*     posterior probability of association from fine-mapping
      - *cs_id:*   ID of 95% credible set (-1 indicates that variant is not in a 95% CS)
      - *beta_posterior:*  posterior expectation of true effect size
      - *sd_posterior:*    posterior standard deviation of true effect size
    - **bgVariants:** path to file of common background variants, with no header and columns `chr`, `start`, `end`, `rsid`. We collated this data from: https://alkesgroup.broadinstitute.org/LDSCORE/baseline_v1.1_hg38_annots/ and the list of background variants can be downloaded here: https://mitra.stanford.edu/kundaje/kmualim/out/Benchmarking/GTEx/all.bg.SNPs.hg38.baseline.v1.1.bed.sorted
    - **nBGVariants:** number of variants in above file
    - **chrSizes:** path to chromosome sizes file for genome build you are using
    - **partition:** path to genome partition .bed file for filtering to distal, noncoding regions with the minimal columns `chr`, `start`, `end`, `category`. The categories "ABC," "AllPeaks," "Other," and "OtherIntron" correspond to distal noncoding regions. To change this, edit the rule `filter_all_variants` in `workflow/rules/preprocessing.smk`.
    - **TSS:** path to .bed file specifying transcription start sites for each gene. The TSS universe file used in our analyses can be downloaded here: https://mitra.stanford.edu/kundaje/kmualim/out/Benchmarking/GTEx/gencode.v29.transcripts.level12.basic.protein_coding.TSS500bp.bed 
    - **distances:** distance stratifications in base pairs of eGene-eVariant separation at which to calculate metrics. Include 30Mb (30000000) to specify all distances. We recommend the using the following distances: [10000, 100000, 2500000, 30000000]
    - **nThresholdSteps:** the number of thresholds to include when generating enrichment-recall curves. We recommend using about 50 steps.
    - **GTExExpression:** path to .gct file from GTEx with transcripts per million (TPM) values for genes in each GTEx tissue. This file is used to filter variants to those linked to genes that are expressed in their respective tissues for the analysis. We downloaded this data from: https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz. The file can also be downloaded here: https://mitra.stanford.edu/kundaje/kmualim/out/Benchmarking/GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct or on Synapse here: https://doi.org/10.7303/syn18134906
    - **GTExGeneUniverse:** path to .bed file specifying genes considered by GTEx. Variants and predictions will be filtered to this set of genes. The gene universe file used in our analysis can be downloaded here: https://mitra.stanford.edu/kundaje/kmualim/out/Benchmarking/GTEx/gencode.v26.genes.collapsed.bed
    - **thresholdTPM:** threshold TPM value to consider an eGene "expressed." We used a threshold of 1 TPM. 
    - **methods:** list of methods to run, with details specified in the methods config table (see below)
    - **GTExTissues:** list of GTEx tissues to consider. We used all that were present in the fine-mapped data: Brain_Cortex, Prostate, Muscle_Skeletal, Artery_Tibial, Skin_Not_Sun_Exposed_Suprapubic, Esophagus_Muscularis, Brain_Cerebellum, Pituitary, Adipose_Subcutaneous, Pancreas, Spleen, Adrenal_Gland, Testis, Lung, Skin_Sun_Exposed_Lower_leg, Heart_Atrial_Appendage, Adipose_Visceral_Omentum, Esophagus_Mucosa, Thyroid, Brain_Nucleus_accumbens_basal_ganglia, Colon_Sigmoid, Breast_Mammary_Tissue, Artery_Aorta, Cells_EBV-transformed_lymphocytes, Heart_Left_Ventricle, Artery_Coronary, Esophagus_Gastroesophageal_Junction, Nerve_Tibial, Liver, Whole_Blood, Colon_Transverse, Stomach
  
2. **Methods config table** - an example file is included at: `config/methods-config-GM12878_predictors.tsv`. The following columns (with a header) are required, where each row represents one predictive method:
   - **method:** predictive method name, to match listed methods in the main config file outlined above
   - **sampleKey:** path to .tsv file with columns and header `biosample`, `predictionFile`, and `GTExTissue`. The `GTExTissue` column should be left EMPTY for biosamples without a match. Multiple biosamples may have the same GTEx tissue match, but and enrichment-recall curve will only be generated for one (the first) biosample match at this time. The prediction files must be separated by biosample and must have a header with minimally the columns `chr`, `start`, `end`, `TargetGene`, `CellType`, `score_col` (see below)
   - **geneUniverse:** path to .bed file specifying which genes are considered by this method. Variants will be filtered to this set of genes for the analysis.
   - **pred_name_long:** a string with the full method name to be used in plots
   - **threshold:** score threshold value to be used to generate boxplots stratified by distance and enrichment heat maps (outputs (1) and (3)). We used the threshold values corresponding to 70% recall from our CRISPR benchmarking pipeline for each model.
   - **score_col:** name of the column in prediction files with the prediction score
   - **color:** hex code specifying what color to use in plots for this method. May be left blank and will be automatically filled.
   - **inverse_predictor:** `TRUE` if high scores correspond to lower prediction confidence (eg. distance to TSS), otherwise `FALSE`
   - **binary:** `TRUE` if this is a binary 0 or 1 predictor, otherwise `FALSE`


