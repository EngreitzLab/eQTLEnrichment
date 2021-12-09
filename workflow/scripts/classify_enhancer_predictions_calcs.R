suppressPackageStartupMessages({
  library(biomaRt)
  library(dplyr)
  library(tidyr)
  library(optparse)
  library(stringr)
  library(gdata)
  library(ggplot2)})

  option_list <- list(
    make_option(c("--variants"), type="character", default=NA, help="tsv with variants"),
    make_option(c("--pred"), type="character", default=NA, help="intersection of variants with enhancer predictions (col 1 = unique ID, col 2 = predicted target gene, col 3 = score"),
    make_option(c("--outDir"), type="character", default=NA, help="output directory for PR plots"),
    make_option(c("--outThresh"), type="character", default=NA, help="outut file for threshold table"),
    make_option(c("--biosample"), type="character", default=NA, help="biosample for enhancer predictions"))
  
  opt = parse_args(OptionParser(option_list=option_list))
  varFile = opt$variants; predictionsIntFile = opt$pred; outDir = opt$outDir; outThresh=opt$outThresh; biosample=opt$biosample
  
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
  
  # plot PR curves
  x=10001
  if (method=="reads_by_dist_to_tss" || method=="reads_by_dist_to_tss_norm"){
    v = data.frame(matrix(0, nrow = x, ncol = 5)) %>% setNames(c("threshold", "prediction.rate.GivenEnhancer",
                                                                    "prediction.rate.total", "prediction.rate.inEnhancer",
                                                                 "nVariants"))
    #v$threshold = seq(from=min(predictionsInt$Score), to=max(predictionsInt$Score), length.out=x)
    v$threshold = seq(from=0, to=0.02, length.out=x)
    for (i in 1:x){
      predTable = dplyr::filter(variants, Score>v$threshold[i])
      v$prediction.rate.GivenEnhancer[i] = nrow(filter(predTable, predictionClass=='inEnhancer-correctGene'))/nrow(filter(predTable,predictionClass!='noOverlap'))
      v$prediction.rate.total[i] = nrow(filter(predTable, predictionClass=='inEnhancer-correctGene'))/nrow(predTable)
      v$prediction.rate.inEnhancer[i] =nrow(filter(predTable,predictionClass!='noOverlap'))/nrow(predTable)
      v$nVariants[i]= nrow(predTable)

    }
  } else if (method=="dist_to_gene" || method=="dist_to_tss"){
    v = data.frame(matrix(0, nrow = x, ncol = 5)) %>% setNames(c("threshold", "prediction.rate.GivenEnhancer",
                                                                    "prediction.rate.total", "prediction.rate.inEnhancer",
                                                                 "nVariants"))
    v$threshold = seq(from=0, to=100000, length.out=x)
    for (i in 1:x){
      predTable = dplyr::filter(variants, Score<v$threshold[i])
      v$prediction.rate.GivenEnhancer[i] = nrow(filter(predTable, predictionClass=='inEnhancer-correctGene'))/nrow(filter(predTable,predictionClass!='noOverlap'))
      v$prediction.rate.total[i] = nrow(filter(predTable, predictionClass=='inEnhancer-correctGene'))/nrow(predTable)
      v$prediction.rate.inEnhancer[i] =nrow(filter(predTable,predictionClass!='noOverlap'))/nrow(predTable)
      v$nVariants[i]= nrow(predTable)
    }

  }
  else {
    v = data.frame(matrix(0, nrow = x, ncol = 5)) %>% setNames(c("threshold", "prediction.rate.GivenEnhancer",
                                                                 "prediction.rate.total", "prediction.rate.inEnhancer",
                                                                 "nVariants"))
    predTable = dplyr::filter(variants, Score==1)
    v$prediction.rate.GivenEnhancer[1] = nrow(filter(predTable, predictionClass=='inEnhancer-correctGene'))/nrow(filter(predTable,predictionClass!='noOverlap'))
    v$prediction.rate.total[1] = nrow(filter(predTable, predictionClass=='inEnhancer-correctGene'))/nrow(predTable)
    v$prediction.rate.inEnhancer[1] =nrow(filter(predTable,predictionClass!='noOverlap'))/nrow(predTable)
    v$nVariants[1]= nrow(predTable)
    
    predTable = dplyr::filter(variants, Score>=0)
    v$prediction.rate.GivenEnhancer[2] = nrow(filter(predTable, predictionClass=='inEnhancer-correctGene'))/nrow(filter(predTable,predictionClass!='noOverlap'))
    v$prediction.rate.total[2] = nrow(filter(predTable, predictionClass=='inEnhancer-correctGene'))/nrow(predTable)
    v$prediction.rate.inEnhancer[2] =nrow(filter(predTable,predictionClass!='noOverlap'))/nrow(predTable)
    v$nVariants[2]= nrow(predTable)
  }

  if (method %in% c("reads_by_dist_to_tss", "reads_by_dist_to_tss_norm", "dist_to_gene", "dist_to_tss")) {
      pr = ggplot(v, aes(x=prediction.rate.total, y=prediction.rate.GivenEnhancer, label=threshold)) +
        geom_point() + geom_text(check_overlap=TRUE, size=5) +
        theme_minimal() + ggtitle(method) + ylim(c(0,1)) + xlim(c(0,1)) +
        theme(text = element_text(size = 20)) +
        xlab('Overlaps predicted enhancer') + ylab('Overlaps predicted enhancer\nlinked to correct gene')

    other = ggplot(v, aes(x=threshold, y=prediction.rate.inEnhancer)) +
      geom_point() +
      theme_minimal() + ggtitle(method) + ylim(c(0,1)) +
      theme(text = element_text(size = 20)) +
      xlab('Threshold') + ylab('Overlaps predicted enhancer')

    #file.name2 = paste0(outDir, "/", "PRPlots/", tissue, ".", method, ".threshold")
    #pdf(file=paste0(file.name2, ".pdf"), width=7,height=5); print(other); dev.off()

   }

  # filter out rows with eGenes connected to un-expressed genes, write prediction table
  fName = file.path(outDir, "thresholdTables", method, paste0(tissue, ".", biosample, ".tsv"))
  write.table(v, file=fName, sep="\t", quote=F, row.names=F, col.names=T)
  write.table(variants, file="", sep="\t", quote=F, row.names=F, col.names=T)
}


write_table()
