## load libraries
library(dplyr)
library(ggplot2)
library(colorspace)

## load data
base_dir = '/oak/stanford/groups/engreitz/Users/sheth/eQTLEnrichment-thresholding/eQTLEnrichment/workflow/'
# methods = c("dist_to_tss", "dist_to_gene", "reads_by_dist_to_tss", 
#             "reads_by_dist_to_tss_norm", "nearest_gene", "nearest_tss", "within_100kb_of_tss",
#             "H3K27ac_reads_by_dist_to_tss", "H3K27ac_reads_by_dist_to_tss_norm",
#             "GraphRegCl", "EpiMap", "ABC_refseq", "Ramil",
#             "GraphReg_ABCScore_PET_CTCF_RamilCorr_numTSSEnhGene",
#             "Top15_shap_scores", "Top20_shap_scores", "Top25_shap_scores")

methods = c("ABC_refseq", "Ramil", "GraphReg_ABCScore_PET_CTCF_RamilCorr_numTSSEnhGene", "Top15_shap_scores", "Top20_shap_scores", "Top25_shap_scores")

# read in tables
dist_to_tss = read.table(file=paste0(base_dir, '2022_0512-dist_to_tss/dist_to_tss/ERCurveTable.tsv'), header=TRUE)
dist_to_gene = read.table(file=paste0(base_dir, '2022_0512-dist_to_gene/dist_to_gene/ERCurveTable.tsv'), header=TRUE)
reads_by_dist = read.table(file=paste0(base_dir, '2022_0512-reads_by_dist/reads_by_dist_to_tss/ERCurveTable.tsv'), header=TRUE)
reads_by_dist_norm = read.table(file=paste0(base_dir, '2022_0512-DHS-norm/reads_by_dist_to_tss_norm/ERCurveTable.tsv'), header=TRUE)
nearest_gene = read.table(file=paste0(base_dir, '2022_0512-nearest_gene/nearest_gene/ERCurveTable.tsv'), header=TRUE)
nearest_tss = read.table(file=paste0(base_dir, '2022_0512-nearest_tss/nearest_tss/ERCurveTable.tsv'), header=TRUE)
within_100kb = read.table(file=paste0(base_dir, '2022_0512-within_100kb/within_100kb_of_tss/ERCurveTable.tsv'), header=TRUE)
H3K27ac_reads_by_dist = read.table(file=paste0(base_dir, '2022_0512-H3K27ac_reads_by_dist/H3K27ac_reads_by_dist_to_tss/ERCurveTable.tsv'), header=TRUE)
H3K27ac_reads_by_dist_norm = read.table(file=paste0(base_dir, '2022_0512-H3K27ac-norm/H3K27ac_reads_by_dist_to_tss_norm/ERCurveTable.tsv'), header=TRUE)
# GraphRegCl.1 = read.table(file=paste0(base_dir, '2022_0629-GraphRegCl/GraphRegCl/ERCurveTable.tsv'), header=TRUE)
# GraphRegCl.2 = read.table(file=paste0(base_dir, '2022_0720-GraphRegCl/GraphRegCl/ERCurveTable.tsv'), header=TRUE)
# GraphRegCl = rbind(GraphRegCl.2, GraphRegCl.1); GraphRegCl = GraphRegCl[order(GraphRegCl$threshold),]
EpiMap = read.table(file=paste0(base_dir, '2022_0629-EpiMap/EpiMap/ERCurveTable.tsv'), header=TRUE)
ABC_refseq.1 = read.table(file=paste0(base_dir, '2022_0705_ABC/ABC_refseq/ERCurveTable.tsv'), header=TRUE)
ABC_refseq.2 = read.table(file=paste0(base_dir, '2022_0804-ABC/ABC_refseq/ERCurveTable.tsv'), header=TRUE)
ABC_refseq = rbind(ABC_refseq.1, ABC_refseq.2)
Ramil = read.table(file=paste0(base_dir, '2022_0629-Ramil/Ramil/ERCurveTable.tsv'), header=TRUE)
GraphReg_ABCScore_PET_CTCF_RamilCorr_numTSSEnhGene.1 = read.table(file=paste0(base_dir, '2022_0726-LogRegClassifier_v1/GraphReg_ABCScore_PET_CTCF_RamilCorr_numTSSEnhGene/ERCurveTable.tsv'), header=TRUE)
Top15_shap_scores.1 = read.table(file=paste0(base_dir, '2022_0726-LogRegClassifier_v1/Top15_shap_scores/ERCurveTable.tsv'), header=TRUE)
Top20_shap_scores.1 = read.table(file=paste0(base_dir, '2022_0726-LogRegClassifier_v1/Top20_shap_scores/ERCurveTable.tsv'), header=TRUE)
Top25_shap_scores.1 = read.table(file=paste0(base_dir, '2022_0726-LogRegClassifier_v1/Top25_shap_scores/ERCurveTable.tsv'), header=TRUE)
GraphReg_ABCScore_PET_CTCF_RamilCorr_numTSSEnhGene.2 = read.table(file=paste0(base_dir, '2022_0802-LogRegClassifier_v1/GraphReg_ABCScore_PET_CTCF_RamilCorr_numTSSEnhGene/ERCurveTable.tsv'), header=TRUE)
Top15_shap_scores.2 = read.table(file=paste0(base_dir, '2022_0802-LogRegClassifier_v1/Top15_shap_scores/ERCurveTable.tsv'), header=TRUE)
Top20_shap_scores.2 = read.table(file=paste0(base_dir, '2022_0802-LogRegClassifier_v1/Top20_shap_scores/ERCurveTable.tsv'), header=TRUE)
Top25_shap_scores.2 = read.table(file=paste0(base_dir, '2022_0802-LogRegClassifier_v1/Top25_shap_scores/ERCurveTable.tsv'), header=TRUE)
GraphReg_ABCScore_PET_CTCF_RamilCorr_numTSSEnhGene = rbind(GraphReg_ABCScore_PET_CTCF_RamilCorr_numTSSEnhGene.1, GraphReg_ABCScore_PET_CTCF_RamilCorr_numTSSEnhGene.2)
Top15_shap_scores = rbind(Top15_shap_scores.1, Top15_shap_scores.2)
Top20_shap_scores = rbind(Top20_shap_scores.1, Top20_shap_scores.2)
Top25_shap_scores = rbind(Top25_shap_scores.1, Top25_shap_scores.2)


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
#GraphRegCl$method="GraphRegCl"
EpiMap$method="EpiMap"
ABC_refseq$method="ABC_refseq"
Ramil$method="Ramil"
GraphReg_ABCScore_PET_CTCF_RamilCorr_numTSSEnhGene$method = "GraphReg_ABCScore_PET_CTCF_RamilCorr_numTSSEnhGene"
Top15_shap_scores$method = "Top15_shap_scores"
Top20_shap_scores$method = "Top20_shap_scores"
Top25_shap_scores$method = "Top25_shap_scores"

# dist_to_tss = dplyr::select(dist_to_tss, prediction.rate.inEnhancer, meanEnrichment, maxEnrichment, method)
# dist_to_gene = dplyr::select(dist_to_gene, prediction.rate.inEnhancer, meanEnrichment, maxEnrichment, method)
# reads_by_dist = dplyr::select(reads_by_dist, prediction.rate.inEnhancer, meanEnrichment, maxEnrichment, method)
# reads_by_dist_norm = dplyr::select(reads_by_dist_norm, prediction.rate.inEnhancer, meanEnrichment, maxEnrichment, method)
# nearest_gene = dplyr::select(nearest_gene, prediction.rate.inEnhancer, meanEnrichment, maxEnrichment, method)
# nearest_tss = dplyr::select(nearest_tss, prediction.rate.inEnhancer, meanEnrichment, maxEnrichment, method)
# within_100kb = dplyr::select(within_100kb, prediction.rate.inEnhancer, meanEnrichment, maxEnrichment, method)
# H3K27ac_reads_by_dist = dplyr::select(H3K27ac_reads_by_dist, prediction.rate.inEnhancer, meanEnrichment, maxEnrichment, method)
# H3K27ac_reads_by_dist_norm = dplyr::select(H3K27ac_reads_by_dist_norm, prediction.rate.inEnhancer, meanEnrichment, maxEnrichment, method)
# GraphRegCl = dplyr::select(GraphRegCl, prediction.rate.inEnhancer, meanEnrichment, maxEnrichment, method)
# EpiMap = dplyr::select(EpiMap, prediction.rate.inEnhancer, meanEnrichment, maxEnrichment, method)
# ABC_refseq = dplyr::select(ABC_refseq, prediction.rate.inEnhancer, meanEnrichment, maxEnrichment, method, enrichment.LCL, recall.LCL)
# Ramil = dplyr::select(Ramil, prediction.rate.inEnhancer, meanEnrichment, maxEnrichment, method, enrichment.LCL, recall.LCL)
# GraphReg_ABCScore_PET_CTCF_RamilCorr_numTSSEnhGene = dplyr::select(GraphReg_ABCScore_PET_CTCF_RamilCorr_numTSSEnhGene, prediction.rate.inEnhancer, meanEnrichment, maxEnrichment, method, enrichment.LCL, recall.LCL)
# Top15_shap_scores = dplyr::select(Top15_shap_scores, prediction.rate.inEnhancer, meanEnrichment, maxEnrichment, method, enrichment.LCL, recall.LCL)
# Top20_shap_scores = dplyr::select(Top20_shap_scores, prediction.rate.inEnhancer, meanEnrichment, maxEnrichment, method, enrichment.LCL, recall.LCL)
# Top25_shap_scores = dplyr::select(Top25_shap_scores, prediction.rate.inEnhancer, meanEnrichment, maxEnrichment, method, enrichment.LCL, recall.LCL)

# df = rbind(dist_to_tss, dist_to_gene, reads_by_dist, reads_by_dist_norm, nearest_gene, nearest_tss, within_100kb, 
#            H3K27ac_reads_by_dist, H3K27ac_reads_by_dist_norm, GraphRegCl, EpiMap, ABC_refseq, Ramil,
#            GraphReg_ABCScore_PET_CTCF_RamilCorr_numTSSEnhGene, Top15_shap_scores, Top20_shap_scores, Top25_shap_scores)
df = rbind(ABC_refseq,Ramil,
           GraphReg_ABCScore_PET_CTCF_RamilCorr_numTSSEnhGene, Top15_shap_scores, Top20_shap_scores, Top25_shap_scores)
df[is.na(df)] = 0

## make color palette
cp = data.frame(method=methods, 
                hex=qualitative_hcl(n=length(methods), palette="Set 2"))
cpList = split(f=cp$method, x=cp$hex)

# prder by threshold
df = df[order(df$threshold),]

# data frame for chosen thresholds
# order of methods: "ABC_refseq", "Ramil", "GraphReg_ABCScore_PET_CTCF_RamilCorr_numTSSEnhGene", "Top15_shap_scores", "Top20_shap_scores", "Top25_shap_scores"
chosen.thresholds = data.frame(method=methods, 
                               threshold=c(0.02, 0, 0.32, 0.32, 0.35, 0.34), 
                               recall.LCL=c(0.17, 0, 0.26, 0.26, 0.258, 0.258), 
                               enrichment.LCL=c(12.3, 100, 4.4, 4.7, 5.8, 5.2))

g=ggplot(data=df, aes(x=recall.LCL, y=enrichment.LCL, color=method)) + 
  geom_path(size=1) +
  geom_point(data=chosen.thresholds, aes(x=recall.LCL, y=enrichment.LCL, color=method), size=4) +
  scale_color_manual(values=cp$hex) + 
  ylim(c(0, 50)) +
  theme_minimal()

pdf(file=("/oak/stanford/groups/engreitz/Users/sheth/eQTLEnrichment-thresholding/eQTLEnrichment/workflow/er_curves_combined_LCLonly_v2.pdf"), width=8, height=5)
   print(g)
dev.off()
g