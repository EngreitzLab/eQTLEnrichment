## load libraries
library(dplyr)
library(ggplot2)
library(colorspace)

## load data
base_dir = '/oak/stanford/groups/engreitz/Users/sheth/eQTLEnrichment-thresholding/eQTLEnrichment/workflow/'

methods_all = c("ABC_refseq", "Ramil", "ABC_Nasser_GM12878_ENCODE", "ABC_Nasser_GM12878_Roadmap", 
            "GraphReg_ABCScore", "GraphReg_ABCScore_RamilWeighted", "GraphReg_ABCScore_PET", "GraphReg_ABCScore_CTCF",
            "GraphReg_ABCScore_PET_CTCF", "GraphReg_ABCScore_PET_CTCF_RamilWeighted", "GraphReg_ABCScore_PET_CTCF_RamilWeighted_sumNearbyEnhancers_10kb",
            "Full_DNase", "Full_DNaseH3K27ac", "Top15_shap_scores", "Top10_shap_scores")

methods = c("ABC_refseq", "Ramil", "EpiMap", "GraphReg", "HiC",
            "Full_v3", "dist_to_tss", "DNase_v3", "Enformer")
# read in tables
dist_to_tss = read.table(file=paste0(base_dir, '2022_0917-dist_to_tss/dist_to_tss/ERCurveTable.tsv'), header=TRUE)
EpiMap = read.table(file=paste0(base_dir, '2022_0830-EpiMap-GM12878/EpiMap_GM12878/ERCurveTable.tsv'), header=TRUE)
ABC_refseq.1 = read.table(file=paste0(base_dir, '2022_0705_ABC/ABC_refseq/ERCurveTable.tsv'), header=TRUE)
ABC_refseq.2 = read.table(file=paste0(base_dir, '2022_0804-ABC/ABC_refseq/ERCurveTable.tsv'), header=TRUE)
ABC_refseq = rbind(ABC_refseq.1, ABC_refseq.2)
Ramil = read.table(file=paste0(base_dir, '2022_1030-EPIraction/Ramil/ERCurveTable.tsv'), header=TRUE)
Full_v3 = read.table(file=paste0(base_dir, '2022_0828-LogRegClassifier_v3_full/full_v3/ERCurveTable.tsv'), header=TRUE)
GraphReg = read.table(file=paste0(base_dir, '2022_0901-new_methods/GraphReg/ERCurveTable.tsv'), header=TRUE)
#HiC = read.table(file=paste0(base_dir, '2022_0916-hic/HiC/ERCurveTable.tsv'), header=TRUE)
ABC_QNORM_hic.1 = read.table(file=paste0(base_dir, '2022_0915-ABC-QNORM-hic/ABC_QNORM-hic/ERCurveTable.tsv'), header=TRUE)
ABC_QNORM_hic.2 = read.table(file=paste0(base_dir, '2022_0916-ABC-QNORM-hic/ABC_QNORM-hic/ERCurveTable.tsv'), header=TRUE)
ABC_QNORM_hic = rbind(ABC_QNORM_hic.1, ABC_QNORM_hic.2)
DNase_v3 = read.table(file=paste0(base_dir, '2022_0908-LogRegClassifier_v3_DNase/v3_DNase/ERCurveTable.tsv'), header=TRUE)
ABC_DNaseOnly = read.table(file=paste0(base_dir, '2022_0917-ABC_DNaseOnly/ABC_DNaseOnly/ERCurveTable.tsv'), header=TRUE)
Enformer = read.table(file=paste0(base_dir, '2022_1030-Enformer/Enformer/ERCurveTable.tsv'), header=TRUE)

# ABC_Nasser_GM12878_ENCODE = read.table(file=paste0(base_dir, '2022_0808-ABC_Nasser/ABC_Nasser_GM12878-ENCODE/ERCurveTable.tsv'), header=TRUE)
# ABC_Nasser_GM12878_Roadmap = read.table(file=paste0(base_dir, '2022_0808-ABC_Nasser/ABC_Nasser_GM12878-Roadmap/ERCurveTable.tsv'), header=TRUE)
# Top10_shap_scores.1 = read.table(file=paste0(base_dir, '2022_0808-LogRegClassifier_v2_full/Top10_shap_scores/ERCurveTable.tsv'), header=TRUE)
# Top15_shap_scores.1 = read.table(file=paste0(base_dir, '2022_0808-LogRegClassifier_v2_full/Top15_shap_scores/ERCurveTable.tsv'), header=TRUE)
# GraphReg_ABCScore.1 = read.table(file=paste0(base_dir, '2022_0808-LogRegClassifier_v2_full/GraphReg+ABCScore/ERCurveTable.tsv'), header=TRUE)
# GraphReg_ABCScore_RamilWeighted.1 = read.table(file=paste0(base_dir, '2022_0808-LogRegClassifier_v2_full/GraphReg+ABCScore+RamilWeighted/ERCurveTable.tsv'), header=TRUE)
# GraphReg_ABCScore_PET.1 = read.table(file=paste0(base_dir, '2022_0808-LogRegClassifier_v2_full/GraphReg+ABCScore+PET/ERCurveTable.tsv'), header=TRUE)
# GraphReg_ABCScore_CTCF.1 = read.table(file=paste0(base_dir, '2022_0808-LogRegClassifier_v2_full/GraphReg+ABCScore+CTCF/ERCurveTable.tsv'), header=TRUE)
# GraphReg_ABCScore_PET_CTCF.1 = read.table(file=paste0(base_dir, '2022_0808-LogRegClassifier_v2_full/GraphReg+ABCScore+PET+CTCF/ERCurveTable.tsv'), header=TRUE)
# GraphReg_ABCScore_PET_CTCF_RamilWeighted.1 = read.table(file=paste0(base_dir, '2022_0808-LogRegClassifier_v2_full/GraphReg+ABCScore+PET+CTCF+RamilWeighted/ERCurveTable.tsv'), header=TRUE)
# GraphReg_ABCScore_PET_CTCF_RamilWeighted_sumNearbyEnhancers_10kb.1 = read.table(file=paste0(base_dir, '2022_0808-LogRegClassifier_v2_full/GraphReg_ABCScore+PET+CTCF+RamilWeighted+sumNearbyEnhancers_10kb/ERCurveTable.tsv'), header=TRUE)
# Full_DNase.1 = read.table(file=paste0(base_dir, '2022_0808-LogRegClassifier_v2_full/Full_DNase/ERCurveTable.tsv'), header=TRUE)
# Full_DNaseH3K27ac.1 = read.table(file=paste0(base_dir, '2022_0808-LogRegClassifier_v2_full/Full_DNaseH3K27ac/ERCurveTable.tsv'), header=TRUE)
# Top10_shap_scores.2 = read.table(file=paste0(base_dir, '2022_0809-LogRegClassifier_v2_full/Top10_shap_scores/ERCurveTable.tsv'), header=TRUE)
# Top15_shap_scores.2 = read.table(file=paste0(base_dir, '2022_0809-LogRegClassifier_v2_full/Top15_shap_scores/ERCurveTable.tsv'), header=TRUE)
# GraphReg_ABCScore.2 = read.table(file=paste0(base_dir, '2022_0809-LogRegClassifier_v2_full/GraphReg+ABCScore/ERCurveTable.tsv'), header=TRUE)
# GraphReg_ABCScore_RamilWeighted.2 = read.table(file=paste0(base_dir, '2022_0809-LogRegClassifier_v2_full/GraphReg+ABCScore+RamilWeighted/ERCurveTable.tsv'), header=TRUE)
# GraphReg_ABCScore_PET.2 = read.table(file=paste0(base_dir, '2022_0809-LogRegClassifier_v2_full/GraphReg+ABCScore+PET/ERCurveTable.tsv'), header=TRUE)
# GraphReg_ABCScore_CTCF.2 = read.table(file=paste0(base_dir, '2022_0809-LogRegClassifier_v2_full/GraphReg+ABCScore+CTCF/ERCurveTable.tsv'), header=TRUE)
# GraphReg_ABCScore_PET_CTCF.2 = read.table(file=paste0(base_dir, '2022_0809-LogRegClassifier_v2_full/GraphReg+ABCScore+PET+CTCF/ERCurveTable.tsv'), header=TRUE)
# GraphReg_ABCScore_PET_CTCF_RamilWeighted.2 = read.table(file=paste0(base_dir, '2022_0809-LogRegClassifier_v2_full/GraphReg+ABCScore+PET+CTCF+RamilWeighted/ERCurveTable.tsv'), header=TRUE)
# GraphReg_ABCScore_PET_CTCF_RamilWeighted_sumNearbyEnhancers_10kb.2 = read.table(file=paste0(base_dir, '2022_0809-LogRegClassifier_v2_full/GraphReg_ABCScore+PET+CTCF+RamilWeighted+sumNearbyEnhancers_10kb/ERCurveTable.tsv'), header=TRUE)
# Full_DNase.2 = read.table(file=paste0(base_dir, '2022_0809-LogRegClassifier_v2_full/Full_DNase/ERCurveTable.tsv'), header=TRUE)
# Full_DNaseH3K27ac.2 = read.table(file=paste0(base_dir, '2022_0809-LogRegClassifier_v2_full/Full_DNaseH3K27ac/ERCurveTable.tsv'), header=TRUE)
# Top10_shap_scores = rbind(Top10_shap_scores.1, Top10_shap_scores.2)
# Top15_shap_scores = rbind(Top15_shap_scores.1, Top15_shap_scores.2)
# GraphReg_ABCScore = rbind(GraphReg_ABCScore.1, GraphReg_ABCScore.2)
# GraphReg_ABCScore_RamilWeighted = rbind(GraphReg_ABCScore_RamilWeighted.1, GraphReg_ABCScore_RamilWeighted.2)
# GraphReg_ABCScore_PET = rbind(GraphReg_ABCScore_PET.1, GraphReg_ABCScore_PET.2)
# GraphReg_ABCScore_CTCF = rbind(GraphReg_ABCScore_CTCF.1, GraphReg_ABCScore_CTCF.2)
# GraphReg_ABCScore_PET_CTCF = rbind(GraphReg_ABCScore_PET_CTCF.1, GraphReg_ABCScore_PET_CTCF.2)
# GraphReg_ABCScore_PET_CTCF_RamilWeighted = rbind(GraphReg_ABCScore_PET_CTCF_RamilWeighted.1, GraphReg_ABCScore_PET_CTCF_RamilWeighted.2)
# GraphReg_ABCScore_PET_CTCF_RamilWeighted_sumNearbyEnhancers_10kb = rbind(GraphReg_ABCScore_PET_CTCF_RamilWeighted_sumNearbyEnhancers_10kb.1, GraphReg_ABCScore_PET_CTCF_RamilWeighted_sumNearbyEnhancers_10kb.2)
# Full_DNase = rbind(Full_DNase.1, Full_DNase.2) 
# Full_DNaseH3K27ac = rbind(Full_DNaseH3K27ac.1, Full_DNaseH3K27ac.2)

dist_to_tss$method="dist_to_tss"
dist_to_tss$plotName="DHS distance to TSS"
# dist_to_gene$method="dist_to_gene"
# reads_by_dist$method="reads_by_dist_to_tss"
# reads_by_dist_norm$method = "reads_by_dist_to_tss_norm"
# ABC$method = "ABC_20210712_hg38ENCODE_RefSeq"
# nearest_gene$method = "nearest_gene"
# nearest_tss$method = "nearest_tss"
# within_100kb$method = "within_100kb_of_tss"
# H3K27ac_reads_by_dist$method = "H3K27ac_reads_by_dist_to_tss"
# H3K27ac_reads_by_dist_norm$method = "H3K27ac_reads_by_dist_to_tss_norm"
# GraphRegCl$method="GraphRegCl"
EpiMap$method="EpiMap"
EpiMap$plotName="EpiMap"
#ABC_refseq$method="ABC_refseq"
#ABC_refseq$plotName="ABC"
ABC_QNORM_hic$method="ABC_QNORM-hic"
ABC_QNORM_hic$plotName="ABC (Intact Hi-C)"
Ramil$method="Ramil"
Ramil$plotName="EPIraction"
Full_v3$method="Full_v3"
Full_v3$plotName="ENCODE-E2G_extended"
GraphReg$method="GraphReg"
GraphReg$plotName="GraphReg_RawScores"
# HiC$method="HiC"
# HiC$plotName="Intact Hi-C interaction freq. (5kb)"
DNase_v3$method="DNase_v3"
DNase_v3$plotName="ENCODE-E2G"
ABC_DNaseOnly$method = "ABC_DNaseOnly"
ABC_DNaseOnly$plotName = "ABC (DNase-only, avg. Intact Hi-C)"
Enformer$method = "Enformer"
Enformer$plotName = "Enformer_Variants"
# GraphReg_ABCScore_PET_CTCF_RamilCorr_numTSSEnhGene$method = "GraphReg_ABCScore_PET_CTCF_RamilCorr_numTSSEnhGene"
# Top10_shap_scores$method = "Top10_shap_scores"
# Top10_shap_scores$plotName = "Log. Regr. (Top 10 features)"
# Top15_shap_scores$method = "Top15_shap_scores"
# Top15_shap_scores$plotName = "Log. Regr. (Top 15 features)"
# Top20_shap_scores$method = "Top20_shap_scores"
# Top25_shap_scores$method = "Top25_shap_scores"
# ABC_Nasser_GM12878_ENCODE$method = "ABC_Nasser_GM12878_ENCODE"
# ABC_Nasser_GM12878_Roadmap$method = "ABC_Nasser_GM12878_Roadmap"
# GraphReg_ABCScore$method = "GraphReg_ABCScore"
# GraphReg_ABCScore_RamilWeighted$method = "GraphReg_ABCScore_RamilWeighted"
# GraphReg_ABCScore_PET$method = "GraphReg_ABCScore_PET"
# GraphReg_ABCScore_CTCF$method = "GraphReg_ABCScore_CTCF"
# GraphReg_ABCScore_PET_CTCF$method = "GraphReg_ABCScore_PET_CTCF"
# GraphReg_ABCScore_PET_CTCF_RamilWeighted$method = "GraphReg_ABCScore_PET_CTCF_RamilWeighted"
# GraphReg_ABCScore_PET_CTCF_RamilWeighted_sumNearbyEnhancers_10kb$method = "GraphReg_ABCScore_PET_CTCF_RamilWeighted_sumNearbyEnhancers_10kb"
# Full_DNase$method = "Full_DNase"
# Full_DNase$plotName = "Log. Regr. (DNase only)"
# Full_DNaseH3K27ac$method = "Full_DNaseH3K27ac"
# Full_DNaseH3K27ac$plotName = "Log. Regr. (DNase + H3K27ac)"

# df_all = rbind(ABC_refseq, Ramil, ABC_Nasser_GM12878_ENCODE, ABC_Nasser_GM12878_Roadmap, 
#            GraphReg_ABCScore, GraphReg_ABCScore_RamilWeighted, GraphReg_ABCScore_PET, GraphReg_ABCScore_CTCF,
#            GraphReg_ABCScore_PET_CTCF, GraphReg_ABCScore_PET_CTCF_RamilWeighted, GraphReg_ABCScore_PET_CTCF_RamilWeighted_sumNearbyEnhancers_10kb,
#            Full_DNase, Full_DNaseH3K27ac, Top15_shap_scores, Top10_shap_scores)

df = rbind(ABC_QNORM_hic,EpiMap, Ramil, Full_v3, dist_to_tss, GraphReg, DNase_v3, ABC_DNaseOnly, Enformer)
df[is.na(df)] = 0

## make color palette
cp = data.frame(method=methods, 
                hex=qualitative_hcl(n=length(methods), palette="Set 2"))
cpList = split(f=cp$method, x=cp$hex)

# prder by threshold
df = df[order(df$threshold),]


# data frame for chosen thresholds
chosen.thresholds = data.frame(method = c("ABC_QNORM_hic","ABC_DNaseOnly", "dist_to_tss", "Ramil", "EpiMap",
                                          "Full_v3", "GraphReg", "DNase_v3",
                                          "nearest_expressed_gene", "nearest_expressed_tss"),
                                plotName = c("ABC (Intact Hi-C)", "ABC (DNase-only, avg. Intact Hi-C)", "DHS distance to TSS", "EPIraction", "EpiMap",
                                             "ENCODE-E2G_extended", "GraphReg_RawScores", "ENCODE-E2G",
                                             "In DHS element & Closest expr. gene","In DHS element & Closest expr. TSS"),
                                threshold = c(0.02, 0.013937, 58811, 0.00032, 0.00452, 0.319, 0,  0.20, 1, 1),
                                recall.LCL = c(.2564, 0.26, 0.145, 0.12,0.11,0.18, 0.49, 0.238, 0.15, 0.15),
                                enrichment.LCL = c(10.89, 5.3, 24,10.2,2.2,27.3, 1.23, 14.5, 33.4, 34.6))
 
 df = dplyr::filter(df, recall.LCL!=0, enrichment.LCL!=0)

g=ggplot(data=df, aes(x=recall.LCL, y=enrichment.LCL, color=plotName)) +
  geom_path(size=1) +
  geom_point(data=chosen.thresholds, aes(x=recall.LCL, y=enrichment.LCL, color=plotName), size=4) +
  #scale_color_manual(values=cp$hex) +
  scale_color_manual(values=c("ABC (Intact Hi-C)"="#4E79A7", "EpiMap"="#E15759", "EPIraction"="darkorchid3",
                              "GraphReg_RawScores" = "seagreen",
                              "ENCODE-E2G_extended" = "#B07AA1", 
                              "DNase reads * distance"="darkseagreen2", "DNase reads * distance (norm.)" = "darkseagreen4",
                              "H3K27ac reads * distance"="violetred1", "H3K27ac reads * distance (norm.)" = "violetred3",
                              "DHS distance to gene"="darkorange4", "DHS distance to TSS" = "darkorange1",
                              "DHS nearest gene"="turquoise1", "DHS nearest TSS"="turquoise4",
                              "DHS within 100kb of TSS"="black",
                              "Intact Hi-C interaction freq. (5kb)"="goldenrod4", 
                              "ENCODE-E2G"="#F7ABE2",
                              "ABC (DNase-only, avg. Intact Hi-C)"="#7DAEE3",
                              "Enformer_Variants" = "deeppink1",
                              "In DHS element & Closest expr. gene"='orangered1',"In DHS element & Closest expr. TSS"='orangered3')) +
  ylab("Enrichment (LCL eQTLs/GM12878 predictions)") + xlab("Recall (LCL eQTLs in GM12878 predictions)") +
  labs(col="Predictor") +
  ylim(c(0, 50)) + xlim(0,0.52) +
  theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8))

ggsave('/oak/stanford/groups/engreitz/Users/sheth/eQTLEnrichment-thresholding/eQTLEnrichment/workflow/main_fig2d.eps', g, width=6, height=4)
ggsave('/oak/stanford/groups/engreitz/Users/sheth/eQTLEnrichment-thresholding/eQTLEnrichment/workflow/main_fig2d.pdf', g, width=6, height=4)



# pdf(file=("/oak/stanford/groups/engreitz/Users/sheth/eQTLEnrichment-thresholding/eQTLEnrichment/workflow/er_curves_combined_mostRecent.pdf"), width=10, height=5)
#    print(g)
# dev.off()
