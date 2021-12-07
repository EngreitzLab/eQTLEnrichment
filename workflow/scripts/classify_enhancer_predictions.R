suppressPackageStartupMessages({
  library(biomaRt)
  library(dplyr)
  library(tidyr)
  library(optparse)
  library(stringr)
  library(gdata)
  library(htmlwidgets)
  library(plotly)})

  option_list <- list(
    make_option(c("--variants"), type="character", default=NA, help="tsv with variants"),
    make_option(c("--pred"), type="character", default=NA, help="intersection of variants with enhancer predictions (col 1 = unique ID, col 2 = predicted target gene, col 3 = score"),
    make_option(c("--outDir"), type="character", default=NA, help="output directory for PR plots"))
  
  opt = parse_args(OptionParser(option_list=option_list))
  varFile = opt$variants; predictionsIntFile = opt$pred; outDir = opt$outDir
  
  variants = read.table(file=varFile,header=FALSE, stringsAsFactor=FALSE); 
  colnames(variants) = c("chr", "start", "end", "unique_id", "GTExTissue", "eGene", "PIP", "TPM")
  #variants$annotations=""
  
  predictionsInt = read.table(file=predictionsIntFile,header=FALSE, stringsAsFactor=FALSE)
  colnames(predictionsInt) = c('unique_id', 'targetGene', 'Score')
  predictionsInt$Score=as.numeric(predictionsInt$Score)
  
  tissue = str_split(predictionsIntFile, '/')[[1]]; tissue = tissue[length(tissue)] %>% str_split(pattern='\\.'); tissue = tissue[[1]][1]
  method = str_split(predictionsIntFile, '/')[[1]]; method = method[length(method)-2]

write_table <- function() {    
  for (i in 1:nrow(variants)) {
    temp = filter(predictionsInt, unique_id==variants$unique_id[i])
    variants$nEnhancers[i] = nrow(temp)
    
    if (nrow(temp)==0 && !(variants$eGene[i] %in% temp$targetGene)) {
      variants$predictionClass[i] = 'noOverlap'
      if  (method=="reads_by_dist_to_tss" || method=="reads_by_dist_to_tss_norm") {
        variants$Score[i] = 10
      } else if (method=="dist_to_gene" || method=="dist_to_tss") {
        variants$Score[i] = 0
      } else {
        variants$Score[i] = 1
      }
      #variants$annotation[i]="if1"
    }
    
    if (nrow(temp)>0 && (variants$eGene[i] %in% temp$targetGene)) {
      variants$predictionClass[i] = 'inEnhancer-correctGene'
      if  (method=="reads_by_dist_to_tss" || method=="reads_by_dist_to_tss_norm") {
        variants$Score[i]=max(temp$Score, na.rm=TRUE)
      }  else if (method=="dist_to_gene" || method=="dist_to_tss") {
        variants$Score[i]=min(temp$Score)
      } else {
        variants$Score[i]=1
      }
      #variants$annotation[i]=paste0(variants$annotation[i], "if2")
    } 
    
    if (nrow(temp)>0 && !(variants$eGene[i] %in% temp$targetGene)) {
      variants$predictionClass[i] = 'inEnhancer-incorrectGene'
      variants$Score[i]=if_else(length(min(temp$Score))==1, min(temp$Score), min(temp$Score)[1])
      #variants$annotation[i]=paste0(variants$annotation[i],"if3")
    }
  }
  
  # # plot PR curves
  # x=1001
  # if (method=="reads_by_dist_to_tss" || method=="reads_by_dist_to_tss_norm"){
  #   v = data.frame(matrix(0, nrow = x, ncol = 5)) %>% setNames(c("threshold", "prediction.rate.GivenEnhancer",
  #                                                                   "prediction.rate.total", "prediction.rate.inEnhancer",
  #                                                                "nVariants"))
  #   #v$threshold = seq(from=min(predictionsInt$Score), to=max(predictionsInt$Score), length.out=x)
  #   v$threshold = seq(from=0, to=1, length.out=x)
  #   for (i in 1:x){
  #     predTable = dplyr::filter(variants, Score>v$threshold[i])
  #     v$prediction.rate.GivenEnhancer[i] = nrow(filter(predTable, predictionClass=='inEnhancer-correctGene'))/nrow(filter(predTable,predictionClass!='noOverlap'))
  #     v$prediction.rate.total[i] = nrow(filter(predTable, predictionClass=='inEnhancer-correctGene'))/nrow(predTable)
  #     v$prediction.rate.inEnhancer[i] =nrow(filter(predTable,predictionClass!='noOverlap'))/nrow(predTable)
  #     v$nVariants[i]= nrow(predTable)
  #     
  #   }
  # }
  # if (method=="dist_to_gene" || method=="dist_to_tss"){
  #   v = data.frame(matrix(0, nrow = x, ncol = 5)) %>% setNames(c("threshold", "prediction.rate.GivenEnhancer",
  #                                                                   "prediction.rate.total", "prediction.rate.inEnhancer",
  #                                                                "nVariants"))
  #   v$threshold = seq(from=0, to=54000, length.out=x)
  #   for (i in 1:x){
  #     predTable = dplyr::filter(variants, Score<v$threshold[i])
  #     v$prediction.rate.GivenEnhancer[i] = nrow(filter(predTable, predictionClass=='inEnhancer-correctGene'))/nrow(filter(predTable,predictionClass!='noOverlap'))
  #     v$prediction.rate.total[i] = nrow(filter(predTable, predictionClass=='inEnhancer-correctGene'))/nrow(predTable)
  #     v$prediction.rate.inEnhancer[i] =nrow(filter(predTable,predictionClass!='noOverlap'))/nrow(predTable)
  #     v$nVariants[i]= nrow(predTable)
  #   }
  #   
  # }
  # 
  # if (method %in% c("reads_by_dist_to_tss", "reads_by_dist_to_tss_norm", "dist_to_gene", "dist_to_tss")) {
  #     pr = ggplot(v, aes(x=prediction.rate.total, y=prediction.rate.GivenEnhancer, label=threshold)) +
  #       geom_point() + geom_label(check_overlap=TRUE) +
  #       theme_minimal() + ggtitle(method) + ylim(c(0,1)) + xlim(c(0,1)) + 
  #       theme(text = element_text(size = 20)) +
  #       xlab('Overlaps predicted enhancer') + ylab('Overlaps predicted enhancer\nlinked to correct gene')
  #   
  #   file.name = paste0(outDir, "/", "PRPlots/", tissue, ".", method)
  #   
  #   other = ggplot(v, aes(x=threshold, y=prediction.rate.total)) +
  #     geom_point() +
  #     theme_minimal() + ggtitle(method) + ylim(c(0,1)) + 
  #     theme(text = element_text(size = 20)) +
  #     xlab('Threshold') + ylab('Overlaps predicted enhancer')
  #   
  #   file.name2 = paste0(outDir, "/", "PRPlots/", tissue, ".", method, ".threshold")
  #   
  #   write.table(v, file=paste0(file.name, ".tsv"), sep="\t", quote=F, row.names=F, col.names=T)
  #   pdf(file=paste0(file.name, ".pdf"), width=7,height=5); print(pr); dev.off()
  #   pdf(file=paste0(file.name2, ".pdf"), width=7,height=5); print(other); dev.off()
  # }

  # filter out rows with eGenes connected to un-expressed genes, write prediction table
  write.table(variants, file="", sep="\t", quote=F, row.names=F, col.names=T)
}


write_table()
