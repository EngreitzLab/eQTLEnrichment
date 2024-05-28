suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(data.table)})

read_GTEx_variants <- function(file_name, method_filter, cs_filter) {
	df = fread(file_name, sep="\t", header=FALSE)
	colnames(df) = c("chr", "start", "end", "variant_id", "varID_hg38", "allele1", "allele2", "cohort", "method", "tissue", "gene_ensembl", "maf", "beta_marginal", "se_marginal", "z", "pip", "cs_id", "beta_posterior", "sd_posterior")
	message("Initial number of variants: ", nrow(df))

	df = dplyr::select(df, chr, start, end, variant_id, method, tissue, gene_ensembl, pip, cs_id) %>%
		dplyr::filter(method==method_filter)
	
	if(cs_filter) {
		df = dplyr::filter(df, cs_id != -1)
	}
	message("After filtering by method and credible set: ", nrow(df))

	df = dplyr::select(df, chr, start, end, variant_id, tissue, gene_ensembl, pip)
	return(df)
}

convert_gene_ids <- function(df, gene_file) {
	genes = fread(gene_file, sep="\t", header=TRUE)
	genes = dplyr::select(genes, name, Ensembl_ID)
	colnames(genes) = c("gene_hgnc", "gene_ensembl_core")
	df$gene_ensembl_core = substr(df$gene_ensembl, 1,15) # get rid of decimals after ensembl IDs

	# join and drop variants with no match
	df = inner_join(df, genes, by="gene_ensembl_core", relationship = "many-to-many")
	df = dplyr::select(df, -c("gene_ensembl", "gene_ensembl_core"))

	message("After filtering to gene universe: ", nrow(df))
	df = distinct(df)

	message("After filtering to unique rows: ", nrow(df))

	return(df)
}

filter_by_expression <- function(df, expression_file, TPM.thresh) {
  col.key = data.frame(tissue=c('Brain_Cortex', 'Prostate', 'Muscle_Skeletal', 'Artery_Tibial', 
  'Skin_Not_Sun_Exposed_Suprapubic', 'Esophagus_Muscularis', 'Brain_Cerebellum', 'Pituitary', 
  'Adipose_Subcutaneous', 'Pancreas', 'Spleen', 'Adrenal_Gland', 'Testis', 'Lung', 
  'Skin_Sun_Exposed_Lower_leg', 'Heart_Atrial_Appendage', 'Adipose_Visceral_Omentum', 
  'Esophagus_Mucosa', 'Thyroid', 'Brain_Nucleus_accumbens_basal_ganglia', 'Colon_Sigmoid', 
  'Breast_Mammary_Tissue', 'Artery_Aorta', 'Cells_EBV-transformed_lymphocytes', 
  'Heart_Left_Ventricle', 'Artery_Coronary', 'Esophagus_Gastroesophageal_Junction', 
  'Nerve_Tibial', 'Liver', 'Whole_Blood', 'Colon_Transverse', 'Stomach'))
  col.key$colNum = c(15, 46, 41, 8, 47, 32, 14, 45, 3, 44, 50, 5, 52, 39, 48, 34, 4, 31, 53, 19, 28, 23, 6, 25, 35, 7, 30, 42, 38, 56, 29, 51)
  
  # read in expression data, select relevant columns, convert to long format with columns for tissue name and TPM
  expr = read.table(expression_file, header=FALSE, skip=3)
  expr = dplyr::select(expr, c(2, col.key$colNum)) %>% setNames(c('gene_hgnc',col.key$tissue)) %>%
  	pivot_longer(cols=-gene_hgnc, names_to='tissue', values_to='TPM')

  # merge variants and expr data, filter to TPM>threshold
  df = left_join(df, expr, by=c('gene_hgnc', 'tissue'), relationship = "many-to-many") %>% dplyr::filter(TPM>TPM.thresh)

	message("After filtering by expression: ", nrow(df))

	df = dplyr::select(df, -TPM)
	return(df)
}

## Process GTEx fine-mapped variants----------------------------------------------------------------------------

# input files
GTEx_variants = "/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/GTExVariants/GTEx_30tissues_hg38.tsv.gz"
GTEx_expression = "/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/GTExVariants/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz"
genes = "/oak/stanford/groups/engreitz/Users/sheth/eQTLEnrichment-integrated/eQTLEnrichment/resources/genome_annotation/CollapsedGeneBounds.hg38.tsv"
TPM_threshold = 1
out_file = paste0("/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/GTExVariants/GTEx_30tissues_hg38.SUSIE.cs.expressedTPM", TPM_threshold, ".tsv..gz")

# read in GTEx file, name columns, filter to fine-mapping method SUSIE and variants in credible set
df = read_GTEx_variants(file_name=GTEx_variants, method_filter="SUSIE", cs_filter=TRUE)

# map ensembl IDs to hgnc
df = convert_gene_ids(df, genes)

# filter to genes expressed in respective tissue
df = filter_by_expression(df, expression_file = GTEx_expression, TPM.thresh = TPM_threshold)

# write to output (columns:  chr, start, end, variant_id, tissue, gene_hgnc, pip)
fwrite(df, file=out_file, col.names=TRUE, sep="\t")
