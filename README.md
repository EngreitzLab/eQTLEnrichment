# eQTL bechmarking pipeline

## Outputs
Running this pipeline will produce:
1. boxplots showing the enrichment of eQTL variants versus common variants, recall (fraction of variants overlapping predicted enhancers), and fraction of variants overlapping a predicted enhnancer linked to the variant's eGene given that the variant overlaps an enhancer. Predicted enhancers are defined as predictions with a score greater than the specified threshold value. For each method, each point represents variants from one GTEx tissue and predicted enhancers from the specified matching biosample. These metrics will be calculated and displayed at the specified distance stratifications of eGene-eVariant separation.
2. enrichment-recall curves for all specified GTEx tissue-biosample pairings, where enrichment and recall (defined here as the fraction of variants overlapping a predicted enhancer linked to at least one of its eGenes) are calculated at the specified number of threshold values yielding even steps across the range of recall for that method.
3. heatmaps showing enrichment at every tissue/biosample pairing for each method

Outputs (1) and (2) are saved as individual pdfs in the /plots directory and can also be viewed in the final eQTL benchmarking report. Note that (1) is only produced if there is more than one GTEx tissue paired with a biosample across the methods. Output (3) can be viewed in the benchmarking report.

## Running the pipeline
To run the pipeline, call each of the rules "first" through "fifth" in consecutive order as shown in the example run-GM12878_predictors.sh file in the /workflow directory. Each rule is dependent on successful completion of the previous ones to successfully run. Note that at this time, rule "fifth," which generates the final benchmarking report, will not run if there is only one tissue/biosample match.

## Config files
The pipeline requires two config files. The required inputs are outlined below.
1. Main config file -- example file: config/config-hg38-GM12878_predictors.yml
The following fields are required:
    - outDir: output directory name
    - codeDir: directory with scripts for running the pipeline (should point to workflow/scripts/)
    - envDir: directory with conda environmentse (should point to workflow/envs)
    - methodsTable: path to methods config file (see below for specifications)
    - GTExVariants: path to file of fine-mapped GTEx variants, with following columns and no header:
          chromosome      chromosome in hg19 coordinates (autosomes only)
          start   start position of variant in hg19 coordinates (0-indexed)
          end     end position of variant in hg19 coordinates (0-indexed)
          variant unique variant identifier in hg19 (chr:pos:ref:alt)
          variant_hg38    unique variant identifier in hg38 (chr:pos:ref:alt)
          allele1 reference allele in hg19 coordinates
          allele2 alternate allele in hg19 coordinates
          cohort  GWAS cohort
          method  fine-mapping method used
          tissue  name of tissue where gene expression was measured
          gene    name of gene of interest for cis-eQTLs
          maf     allele frequency of the minor allele in cohort
          beta_marginal   marginal association effect size from linear mixed model
          se_marginal     standard error on marginal association effect size from linear mixed model
          z       test statistic for marginal association
          pip     posterior probability of association from fine-mapping
          cs_id   ID of 95% credible set (-1 indicates that variant is not in a 95% CS)
          beta_posterior  posterior expectation of true effect size
          sd_posterior    posterior standard deviation of true effect size
    - bgVariants: path to file of common background variants, with no header and columns chr, start, end, rsid. We downloaded this data from: https://alkesgroup.broadinstitute.org/LDSCORE/baseline_v1.1_hg38_annots/
    - nBGVariants: number of variants in above file
    - chrSizes: path to chromosome sizes file for genome build you are using
    - partition: path to genome partition .bed file for filtering to distal, noncoding regions with the minimal columns chr, start, end, category. The categories "ABC," "AllPeaks," "Other," and "OtherIntron" correspond to distal noncoding regions. To change this, edit the rule filter_all_variants in workflow/rules/preprocessing.smk.
    - TSS: path to .bed file specifying transcription start sites for each gene
    - distances: distance stratifications in base pairs of eGene-eVariant separation to calculate metrics at. Include 30Mb (30000000) to specify all distances. We recommend the using the following distances: [10000, 100000, 2500000, 30000000]
    - nThresholdSteps: the number of thresholds to include when generating enrichment-recall curves. We recommend using about 50 steps.
    - GTExExpression: path to .gct file from GTEx with transcipt per million (TPM) values for genes in each GTEx tissue. This file is used to filter variants to those linked to genes that are expressed in their respective tissues for the analysis. We downloaded this data from: https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz
    - GTExGeneUniverse: path to .bed file specifying genes considered by GTEx. Variants and predictions will be filtered to this set of genes.
    - thresholdTPM: threshold TPM value to consider an eGene "expressed." We used a threshold of 1. 
    - methods: list of methods to run, with details specified in the methods config table (see below)
    - GTExTissues: list of GTEx tissues to consider. We used all that were present in the fine-mapped data: [Brain_Cortex, Prostate, Muscle_Skeletal, Artery_Tibial, Skin_Not_Sun_Exposed_Suprapubic, Esophagus_Muscularis, Brain_Cerebellum, Pituitary, Adipose_Subcutaneous, Pancreas, Spleen, Adrenal_Gland, Testis, Lung, Skin_Sun_Exposed_Lower_leg, Heart_Atrial_Appendage, Adipose_Visceral_Omentum, Esophagus_Mucosa, Thyroid, Brain_Nucleus_accumbens_basal_ganglia, Colon_Sigmoid, Breast_Mammary_Tissue, Artery_Aorta, Cells_EBV-transformed_lymphocytes, Heart_Left_Ventricle, Artery_Coronary, Esophagus_Gastroesophageal_Junction, Nerve_Tibial, Liver, Whole_Blood, Colon_Transverse, Stomach]
  
2. Methods config table -- example file: config/methods-config-GM12878_predictors.tsv
The following columns (with header) are required, where each row represents one predictive method:
   - method: predictive method name, to match listed methods in .yml above
   - sampleKey: path to .tsv file with columns and header "biosample", "predictionFile", and "GTExTissue." The GTExTissue column should be left EMPTY for biosamples without a match. Multiple biosamples may have the same GTEx tissue match, but and enrichment-recall curve will only be generated for one (the first) biosample match at this time. The prediction file must have a header with minimally the columns chr, start, end, TargetGene, CellType, score_col (see below)
   - geneUniverse: path to .bed file specifying which genes are considered by this method. Variants will be filtered to this set of genes for the analysis.
   - pred_name_long: a string with the full method name to be used in plots
   - threshold: score threshold value to be used to generate boxplots stratified by distance and enrichment heat maps (outputs (1) and (3)). We used the threhsold values corresponding to 70% recall from our CRISPR benchmarking pipeline for each model.
   - score_col: name of the column in prediction files with the prediction score
   - color: hex code specifying what color to use in plots for this method. May be left blank and will be automatically filled.
   - inverse_predictor: TRUE if high scores correspond to lower prediction confidence (eg. distance to TSS), otherwise FALSE
   - binary: TRUE if this is a binary 0 or 1 predictor, otherwise FALSE


