## load libraries
library(dplyr)
library(ggplot2)
library(colorspace)

## load data
base_dir = '/oak/stanford/groups/engreitz/Users/sheth/eQTLEnrichment-thresholding/eQTLEnrichment/workflow/'
methods = c("ABC_20210712_hg38ENCODE_RefSeq", "dist_to_tss", "dist_to_gene", "reads_by_dist_to_tss", "
            reads_by_dist_to_tss_norm", "nearest_gene", "nearest_tss", "within_100kb_of_tss")
# read in tables
dist_to_tss = read.table(file=paste0(base_dir, '2022_0118_test2-dist_to_tss/dist_to_tss/ERCurveTable.tsv'), header=TRUE)
dist_to_gene = read.table(file=paste0(base_dir, '2022_0118_test2-dist_to_gene/dist_to_gene/ERCurveTable.tsv'), header=TRUE)
reads_by_dist = read.table(file=paste0(base_dir, '2022_0118_test2-reads_by_dist/reads_by_dist_to_tss/ERCurveTable.tsv'), header=TRUE)
reads_by_dist_norm = read.table(file=paste0(base_dir, '2022_0109_test2-norm/reads_by_dist_to_tss_norm/ERCurveTable.tsv'), header=TRUE)
ABC = read.table(file=paste0(base_dir, '2022_0109_test2-ABC/ABC_20210712_hg38ENCODE_RefSeq/ERCurveTable.tsv'), header=TRUE)
nearest_gene = read.table(file=paste0(base_dir, '2022_0109_test2-nearest_gene/nearest_gene/ERCurveTable.tsv'), header=TRUE)
nearest_tss = read.table(file=paste0(base_dir, '2022_0109_test2-nearest_tss/nearest_tss/ERCurveTable.tsv'), header=TRUE)
within_100kb = read.table(file=paste0(base_dir, '2022_0109_test2-within_100kb/within_100kb_of_tss/ERCurveTable.tsv'), header=TRUE)

dist_to_tss$method="dist_to_tss"
dist_to_gene$method="dist_to_gene"
reads_by_dist$method="reads_by_dist_to_tss"
reads_by_dist_norm$method = "reads_by_dist_to_tss_norm"
ABC$method = "ABC_20210712_hg38ENCODE_RefSeq"
nearest_gene$method = "nearest_gene"
nearest_tss$method = "nearest_tss"
within_100kb$method = "within_100kb_of_tss"

df = rbind(dist_to_tss, dist_to_gene, reads_by_dist, reads_by_dist_norm, ABC, nearest_gene, nearest_tss, within_100kb)

## make color palette
cp = data.frame(method=methods, 
                hex=qualitative_hcl(n=length(methods), palette="Set 2"))
cpList = split(f=cp$method, x=cp$hex)

g=ggplot(data=df, aes(x=prediction.rate.inEnhancer, y=meanEnrichment, color=method)) + 
  geom_path(size=1.25) +
  geom_point(size=1.25) +
  scale_color_manual(values=cp$hex) + 
  ylim(c(0, 30)) +
  theme_minimal()

pdf(file=paste0(base_dir, '2022_0109_test2/er_curves_combined_v2.pdf'), width=7, height=5)
   print(g)
dev.off()
g