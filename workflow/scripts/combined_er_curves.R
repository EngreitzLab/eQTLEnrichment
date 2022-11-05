## load libraries
library(dplyr)
library(ggplot2)
library(colorspace)

## load data
base_dir = '/oak/stanford/groups/engreitz/Users/sheth/eQTLEnrichment-thresholding/eQTLEnrichment/workflow/'
methods = c("dist_to_tss", "dist_to_gene", "reads_by_dist_to_tss", 
            "reads_by_dist_to_tss_norm", "nearest_gene", "nearest_tss", "within_100kb_of_tss",
            "H3K27ac_reads_by_dist_to_tss", "H3K27ac_reads_by_dist_to_tss_norm",
            "GraphRegCl", "EpiMap", "ABC_refseq", "Ramil")
# read in tables
# dist_to_tss = read.table(file=paste0(base_dir, '2022_0118_test2-dist_to_tss/dist_to_tss/ERCurveTable.tsv'), header=TRUE)
# dist_to_gene = read.table(file=paste0(base_dir, '2022_0118_test2-dist_to_gene/dist_to_gene/ERCurveTable.tsv'), header=TRUE)
# reads_by_dist.1 = read.table(file=paste0(base_dir, '2022_0118_test2-reads_by_dist/reads_by_dist_to_tss/ERCurveTable.2.tsv'), header=TRUE)
# reads_by_dist.2 = read.table(file=paste0(base_dir, '2022_0118_test2-reads_by_dist/reads_by_dist_to_tss/ERCurveTable.10.tsv'), header=TRUE)
# reads_by_dist = read.table(file=paste0(base_dir, '2022_0118_test2-reads_by_dist/reads_by_dist_to_tss/ERCurveTable.tsv'), header=TRUE)
# 
# #reads_by_dist = rbind(reads_by_dist.1, reads_by_dist.2, reads_by_dist.3)
# #reads_by_dist = reads_by_dist[order(as.numeric(reads_by_dist$threshold)), ]
# 
# reads_by_dist_norm = read.table(file=paste0(base_dir, '2022_0109_test2-norm/reads_by_dist_to_tss_norm/ERCurveTable.tsv'), header=TRUE)
# ABC = read.table(file=paste0(base_dir, '2022_0109_test2-ABC/ABC_20210712_hg38ENCODE_RefSeq/ERCurveTable.tsv'), header=TRUE)
# nearest_gene = read.table(file=paste0(base_dir, '2022_0109_test2-nearest_gene/nearest_gene/ERCurveTable.tsv'), header=TRUE)
# nearest_tss = read.table(file=paste0(base_dir, '2022_0109_test2-nearest_tss/nearest_tss/ERCurveTable.tsv'), header=TRUE)
# within_100kb = read.table(file=paste0(base_dir, '2022_0109_test2-within_100kb/within_100kb_of_tss/ERCurveTable.tsv'), header=TRUE)

# new path names
dist_to_tss = read.table(file=paste0(base_dir, '2022_0512-dist_to_tss/dist_to_tss/ERCurveTable.tsv'), header=TRUE)
dist_to_gene = read.table(file=paste0(base_dir, '2022_0512-dist_to_gene/dist_to_gene/ERCurveTable.tsv'), header=TRUE)
reads_by_dist = read.table(file=paste0(base_dir, '2022_0512-reads_by_dist/reads_by_dist_to_tss/ERCurveTable.tsv'), header=TRUE)
reads_by_dist_norm = read.table(file=paste0(base_dir, '2022_0512-DHS-norm/reads_by_dist_to_tss_norm/ERCurveTable.tsv'), header=TRUE)
nearest_gene = read.table(file=paste0(base_dir, '2022_0512-nearest_gene/nearest_gene/ERCurveTable.tsv'), header=TRUE)
nearest_tss = read.table(file=paste0(base_dir, '2022_0512-nearest_tss/nearest_tss/ERCurveTable.tsv'), header=TRUE)
within_100kb = read.table(file=paste0(base_dir, '2022_0512-within_100kb/within_100kb_of_tss/ERCurveTable.tsv'), header=TRUE)
H3K27ac_reads_by_dist = read.table(file=paste0(base_dir, '2022_0512-H3K27ac_reads_by_dist/H3K27ac_reads_by_dist_to_tss/ERCurveTable.tsv'), header=TRUE)
H3K27ac_reads_by_dist_norm = read.table(file=paste0(base_dir, '2022_0512-H3K27ac-norm/H3K27ac_reads_by_dist_to_tss_norm/ERCurveTable.tsv'), header=TRUE)
GraphRegCl.1 = read.table(file=paste0(base_dir, '2022_0629-GraphRegCl/GraphRegCl/ERCurveTable.tsv'), header=TRUE)
GraphRegCl.2 = read.table(file=paste0(base_dir, '2022_0720-GraphRegCl/GraphRegCl/ERCurveTable.tsv'), header=TRUE)
GraphRegCl = rbind(GraphRegCl.2, GraphRegCl.1); GraphRegCl = GraphRegCl[order(GraphRegCl$threshold),]
EpiMap = read.table(file=paste0(base_dir, '2022_0629-EpiMap/EpiMap/ERCurveTable.tsv'), header=TRUE)
ABC_refseq = read.table(file=paste0(base_dir, '2022_0705_ABC/ABC_refseq/ERCurveTable.tsv'), header=TRUE)
Ramil = read.table(file=paste0(base_dir, '2022_0629_Ramil/Ramil/ERCurveTable.tsv'), header=TRUE)

dist_to_tss$method="dist_to_tss"
dist_to_gene$method="dist_to_gene"
reads_by_dist$method="reads_by_dist_to_tss"
reads_by_dist_norm$method = "reads_by_dist_to_tss_norm"
#ABC$method = "ABC_20210712_hg38ENCODE_RefSeq"
nearest_gene$method = "nearest_gene"
nearest_tss$method = "nearest_tss"
within_100kb$method = "within_100kb_of_tss"
H3K27ac_reads_by_dist$method = "H3K27ac_reads_by_dist_to_tss"
H3K27ac_reads_by_dist_norm$method = "H3K27ac_reads_by_dist_to_tss_norm"
GraphRegCl$method="GraphRegCl"
EpiMap$method="EpiMap"
ABC_refseq$method="ABC_refseq"
Ramil$method="Ramil"

dist_to_tss = dplyr::select(dist_to_tss, prediction.rate.inEnhancer, meanEnrichment, maxEnrichment, method)
dist_to_gene = dplyr::select(dist_to_gene, prediction.rate.inEnhancer, meanEnrichment, maxEnrichment, method)
reads_by_dist = dplyr::select(reads_by_dist, prediction.rate.inEnhancer, meanEnrichment, maxEnrichment, method)
reads_by_dist_norm = dplyr::select(reads_by_dist_norm, prediction.rate.inEnhancer, meanEnrichment, maxEnrichment, method)
#ABC = dplyr::select(ABC, prediction.rate.inEnhancer, meanEnrichment, maxEnrichment, method)
nearest_gene = dplyr::select(nearest_gene, prediction.rate.inEnhancer, meanEnrichment, maxEnrichment, method)
nearest_tss = dplyr::select(nearest_tss, prediction.rate.inEnhancer, meanEnrichment, maxEnrichment, method)
within_100kb = dplyr::select(within_100kb, prediction.rate.inEnhancer, meanEnrichment, maxEnrichment, method)
H3K27ac_reads_by_dist = dplyr::select(H3K27ac_reads_by_dist, prediction.rate.inEnhancer, meanEnrichment, maxEnrichment, method)
H3K27ac_reads_by_dist_norm = dplyr::select(H3K27ac_reads_by_dist_norm, prediction.rate.inEnhancer, meanEnrichment, maxEnrichment, method)
GraphRegCl = dplyr::select(GraphRegCl, prediction.rate.inEnhancer, meanEnrichment, maxEnrichment, method)
EpiMap = dplyr::select(EpiMap, prediction.rate.inEnhancer, meanEnrichment, maxEnrichment, method)
ABC_refseq = dplyr::select(ABC_refseq, prediction.rate.inEnhancer, meanEnrichment, maxEnrichment, method)
Ramil = dplyr::select(Ramil, prediction.rate.inEnhancer, meanEnrichment, maxEnrichment, method)

df = rbind(dist_to_tss, dist_to_gene, reads_by_dist, reads_by_dist_norm, nearest_gene, nearest_tss, within_100kb, 
           H3K27ac_reads_by_dist, H3K27ac_reads_by_dist_norm, GraphRegCl, EpiMap, ABC_refseq, Ramil)
df[is.na(df)] = 0

## make color palette
cp = data.frame(method=methods, 
                hex=qualitative_hcl(n=length(methods), palette="Set 2"))
cpList = split(f=cp$method, x=cp$hex)

g=ggplot(data=df, aes(x=prediction.rate.inEnhancer, y=meanEnrichment, color=method)) + 
  geom_path(size=1.25) +
  geom_point(size=1.25) +
  scale_color_manual(values=cp$hex) + 
  ylim(c(0, 35)) +
  theme_minimal()

pdf(file=("/oak/stanford/groups/engreitz/Users/sheth/eQTLEnrichment-thresholding/eQTLEnrichment/workflow/er_curves_combined_v3.pdf"), width=7, height=5)
   print(g)
dev.off()
g