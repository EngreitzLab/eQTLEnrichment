## load libraries
library(dplyr)
library(ggplot2)
library(colorspace)

## load data
base_dir = '/oak/stanford/groups/engreitz/Users/sheth/eQTLEnrichment-thresholding/eQTLEnrichment/workflow/'

#methods = c("dist_to_tss", "dist_to_gene", "reads_by_dist", "reads_by_dist_norm", "nearest_gene", "nearest_tss", "within_100kb_of_tss", "H3K27ac_reads_by_dist", "H3K27ac_reads_by_dist_norm", "ABC_refseq")

# read in tables
#dist_to_tss = read.table(file=paste0(base_dir, '2022_0803-baseline_predictors/dist_to_tss/ERCurveTable.tsv'), header=TRUE)
#dist_to_gene = read.table(file=paste0(base_dir, '2022_0803-baseline_predictors/dist_to_gene/ERCurveTable.tsv'), header=TRUE)
dist_to_tss = read.table(file=paste0(base_dir, '2022_0917-dist_to_tss/dist_to_tss/ERCurveTable.tsv'), header=TRUE)
dist_to_gene = read.table(file=paste0(base_dir, '2022_0917-dist_to_gene/dist_to_gene/ERCurveTable.tsv'), header=TRUE)
reads_by_dist = read.table(file=paste0(base_dir, '2022_0828-reads_by_dist/reads_by_dist_to_tss/ERCurveTable.tsv'), header=TRUE)
reads_by_dist_norm = read.table(file=paste0(base_dir, '2022_0803-baseline_predictors/reads_by_dist_to_tss_norm/ERCurveTable.tsv'), header=TRUE)
nearest_gene = read.table(file=paste0(base_dir, '2022_0915-nearest_gene/nearest_gene/ERCurveTable.tsv'), header=TRUE)
nearest_tss = read.table(file=paste0(base_dir, '2022_0915-nearest_tss/nearest_tss/ERCurveTable.tsv'), header=TRUE)
within_100kb = read.table(file=paste0(base_dir, '2022_0915-within_100kb/within_100kb_of_tss/ERCurveTable.tsv'), header=TRUE)
H3K27ac_reads_by_dist = read.table(file=paste0(base_dir, '2022_0828-H3K27ac_reads_by_dist/H3K27ac_reads_by_dist_to_tss/ERCurveTable.tsv'), header=TRUE)
H3K27ac_reads_by_dist_norm = read.table(file=paste0(base_dir, '2022_0803-baseline_predictors/H3K27ac_reads_by_dist_to_tss_norm/ERCurveTable.tsv'), header=TRUE)
ABC_QNORM_hic.1 = read.table(file=paste0(base_dir, '2022_0915-ABC-QNORM-hic/ABC_QNORM-hic/ERCurveTable.tsv'), header=TRUE)
ABC_QNORM_hic.2 = read.table(file=paste0(base_dir, '2022_0916-ABC-QNORM-hic/ABC_QNORM-hic/ERCurveTable.tsv'), header=TRUE)
ABC_QNORM_hic = rbind(ABC_QNORM_hic.1, ABC_QNORM_hic.2)


dist_to_tss$method="dist_to_tss"
dist_to_tss$plotName="DHS distance to TSS"
dist_to_gene$method="dist_to_gene"
dist_to_gene$plotName="DHS distance to gene"
reads_by_dist$method="reads_by_dist_to_tss"
reads_by_dist$plotName="DNase reads * distance"
reads_by_dist_norm$method = "reads_by_dist_to_tss_norm"
reads_by_dist_norm$plotName = "DNase reads * distance (norm.)"
nearest_gene$method = "nearest_gene"
nearest_gene$plotName = "DHS nearest gene"
nearest_tss$method = "nearest_tss"
nearest_tss$plotName = "DHS nearest TSS"
within_100kb$method = "within_100kb_of_tss"
within_100kb$plotName = "DHS within 100kb of TSS"
H3K27ac_reads_by_dist$method = "H3K27ac_reads_by_dist_to_tss"
H3K27ac_reads_by_dist$plotName = "H3K27ac reads * distance"
H3K27ac_reads_by_dist_norm$method = "H3K27ac_reads_by_dist_to_tss_norm"
H3K27ac_reads_by_dist_norm$plotName = "H3K27ac reads * distance (norm.)"
ABC_QNORM_hic$method="ABC_QNORM_hic"
ABC_QNORM_hic$plotName="ABC (Intact Hi-C)"



df = rbind(dist_to_tss, dist_to_gene, reads_by_dist, reads_by_dist_norm, H3K27ac_reads_by_dist, H3K27ac_reads_by_dist_norm, ABC_QNORM_hic)
df[is.na(df)] = 0

## make color palette
# cp = data.frame(method=methods, 
#                 hex=qualitative_hcl(n=length(methods), palette="Set 2"))
# cpList = split(f=cp$method, x=cp$hex)

#order by threshold
df = df[order(df$threshold),]

# data frame for chosen thresholds

# add points for all methods
 chosen.thresholds = data.frame(method=c("ABC_QNORM_hic", "nearest_tss", "nearest_gene", "within_100kb_of_tss", "nearest_expressed_gene", "nearest_expressed_tss", "dist_to_tss", "dist_to_gene", "H3K27ac_reads_by_dist_to_tss", "H3K27ac_reads_by_dist_to_tss_norm", "reads_by_dist_to_tss", "reads_by_dist_to_tss_norm"), 
                                plotName=c("ABC (Intact Hi-C)", "DHS nearest TSS", "DHS nearest gene", "DHS within 100kb of TSS", "DHS nearest expr. gene", "DHS nearest expr. TSS", "DHS distance to TSS", "DHS distance to gene", "H3K27ac reads * distance", "H3K27ac reads * distance (norm.)", "DNase reads * distance", "DNase reads * distance (norm.)"),
                                threshold=c(0.02,1,1,1,1,1, 53966, 44840, .00022, .00082, .000314, .000739), 
                                recall.LCL=c(.2564, 0.146, 0.158, 0.154, 0.15, 0.15, .145, .15, .156, .156, .156, .156), 
                                enrichment.LCL=c(10.89, 36.11, 34.52, 13.68, 33.4, 34.6, 24, 20, 1.6, .64, 1.54, .65))
 
 df = dplyr::filter(df, recall.LCL!=0, enrichment.LCL!=0)

g=ggplot(data=df, aes(x=recall.LCL, y=enrichment.LCL, color=plotName)) +
  geom_path(size=1) +
  geom_point(data=chosen.thresholds, aes(x=recall.LCL, y=enrichment.LCL, color=plotName), size=4) +
  #scale_color_manual(values=cp$hex) +
  scale_color_manual(values=c("ABC (Intact Hi-C)"="#4E79A7", "EpiMap"="#E15759", "Ramil"="darkorchid3",
                     "Enformer (CRISPR)"="deeppink1", "GraphReg (scores)" = "seagreen",
                     "CIA" = "cyan3", "Log. Regr. (Full model)" = "#B07AA1", 
                     "DNase reads * distance"="darkseagreen2", "DNase reads * distance (norm.)" = "darkseagreen4",
                     "H3K27ac reads * distance"="violetred1", "H3K27ac reads * distance (norm.)" = "violetred3",
                     "DHS distance to gene"="darkorange4", "DHS distance to TSS" = "darkorange1",
                     "DHS nearest gene"="turquoise1", "DHS nearest TSS"="turquoise4",
                     "DHS within 100kb of TSS"="black",
                     "DHS nearest expr. gene"="orangered1", "DHS nearest expr. TSS"="orangered3")) +
  ylab("Enrichment (LCL eQTLs/GM12878 predictions)") + xlab("Recall (LCL eQTLs in GM12878 predictions)") +
  labs(col="Predictor") +
  ylim(c(0, 50)) +
  theme_minimal()

pdf(file=("/oak/stanford/groups/engreitz/Users/sheth/eQTLEnrichment-thresholding/eQTLEnrichment/workflow/er_curves_combined_LCLonly_bp_v2.pdf"), width=8, height=5)
   print(g)
dev.off()
g