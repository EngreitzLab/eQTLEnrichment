suppressPackageStartupMessages({
  library(biomaRt)
  library(dplyr)
  library(tidyr)
  library(optparse)
  library(stringr)
  library(gdata)
  library(ggplot2)})

  ## INPUTS
  # read in arguments
  option_list <- list(
    make_option(c("--variants"), type="character", default=NA, help="tsv with GTEx variants"),
    make_option(c("--pred"), type="character", default=NA, help="intersection of variants with enhancer predictions (col 1 = unique ID, col 2 = predicted target gene, col 3 = score"),
    make_option(c("--outDir"), type="character", default=NA, help="output directory for PR plots"),
    make_option(c("--outThresh"), type="character", default=NA, help="outut file for threshold table"),
    make_option(c("--biosample"), type="character", default=NA, help="biosample for enhancer predictions"),
    make_option(c("--span"), type="character", default=NA, help="threshold span"))
  
  # make variables for arguments
  opt = parse_args(OptionParser(option_list=option_list))
  varFile = opt$variants; predictionsIntFile = opt$pred; outDir = opt$outDir; outThresh=opt$outThresh; biosample=opt$biosample
  span = opt$span %>% strsplit(" ") %>% unlist()
  
  # read in variants file, set column names
  variants = read.table(file=varFile,header=FALSE, stringsAsFactor=FALSE); 
  colnames(variants) = c("chr", "start", "end", "unique_id", "GTExTissue", "eGene", "PIP", "TPM")

  # read in predictions int file, set column names, set score to numeric
  predictionsInt = read.table(file=predictionsIntFile,header=FALSE, stringsAsFactor=FALSE)
  colnames(predictionsInt) = c('unique_id', 'targetGene', 'Score')
  predictionsInt$Score=as.numeric(predictionsInt$Score)
  
  # extract tissue  name and method name from predictions int file name
  tissue = str_split(predictionsIntFile, '/')[[1]]; tissue = tissue[length(tissue)] %>% str_split(pattern='\\.'); tissue = tissue[[1]][1]
  method = str_split(predictionsIntFile, '/')[[1]]; method = method[length(method)-2]

  ## GENERATE PREDICTION TABLE (for unthresholded predictors)
  # Iterate through variants
  for (i in 1:nrow(variants)) {
    temp = filter(predictionsInt, unique_id==variants$unique_id[i]) # filter predictions to given variant
    variants$nEnhancers[i] = nrow(temp) # record number of predictions overlapping this variant
    
    # if this variant does not overlap a prediction (and thus does not have a correct target gene match)
    # set score for prediction to make sure it's always included as incorrectly predicted 
    if (nrow(temp)==0 && !(variants$eGene[i] %in% temp$targetGene)) {
      variants$predictionClass[i] = 'noOverlap'
      if  (method=="reads_by_dist_to_tss" || method=="reads_by_dist_to_tss_norm") {
        variants$Score[i] = 10
      } else if (method=="dist_to_gene" || method=="dist_to_tss") {
        variants$Score[i] = 0
      } else {
        variants$Score[i] = 5
      }
    }
    
    # if the variant does overlap a prediction and is linked to the correct target gene in any of them
    # set score to the "best" of the predictions overlapping the variant
    if (nrow(temp)>0 && (variants$eGene[i] %in% temp$targetGene)) {
      variants$predictionClass[i] = 'inEnhancer-correctGene'
      if (method=="dist_to_gene" || method=="dist_to_tss") {
        variants$Score[i]=min(temp$Score)
      } else {
        variants$Score[i]=max(temp$Score)
      }
    } 
    
    # if the variant overlaps predictions but is not linked to the correct target gene
    # set corresponding score to the "best" score of predictions overlapping variant, as above
    if (nrow(temp)>0 && !(variants$eGene[i] %in% temp$targetGene)) {
      variants$predictionClass[i] = 'inEnhancer-incorrectGene'
      if (method=="dist_to_gene" || method=="dist_to_tss") {
        variants$Score[i]=min(temp$Score)
      } else {
        variants$Score[i]=max(temp$Score)
      }
    }
  }
  write.table(variants, file="", sep="\t", quote=F, row.names=F, col.names=T)

  ## GENERATE THRESHOLD TABLE
  # initialize data frame with columns for threshold table; one row per threshold in the threshold span
  x = length(span)
  v = data.frame(matrix(0, nrow = x, ncol = 5)) %>% setNames(c("threshold", "prediction.rate.GivenEnhancer",
                                                               "prediction.rate.total", "prediction.rate.inEnhancer",
                                                               "nVariants"))
  v$threshold=as.numeric(span) # set threshold values to numerics
  
  # for predictors where high score is "bad" (distance)
  if (method=="dist_to_gene" || method=="dist_to_tss"){
    # iterate across thresholds
    for (i in 1:x){
      predTable = dplyr::filter(variants, Score<v$threshold[i])
      v$prediction.rate.GivenEnhancer[i] = nrow(filter(predTable, predictionClass=='inEnhancer-correctGene'))/nrow(filter(predTable,predictionClass!='noOverlap'))
      v$prediction.rate.total[i] = nrow(filter(predTable, predictionClass=='inEnhancer-correctGene'))/nrow(predTable)
      v$prediction.rate.inEnhancer[i] =nrow(filter(predTable,predictionClass!='noOverlap'))/nrow(predTable)
      v$nVariants[i]= nrow(predTable)
    }
  # for the rest of the predictors, where low score is bad, do the same
  } else { 
    for (i in 1:x){
      predTable = dplyr::filter(variants, Score>=v$threshold[i])
      v$prediction.rate.GivenEnhancer[i] = nrow(filter(predTable, predictionClass=='inEnhancer-correctGene'))/nrow(filter(predTable,predictionClass!='noOverlap'))
      v$prediction.rate.total[i] = nrow(filter(predTable, predictionClass=='inEnhancer-correctGene'))/nrow(predTable)
      v$prediction.rate.inEnhancer[i] =nrow(filter(predTable,predictionClass!='noOverlap'))/nrow(predTable)
      v$nVariants[i]= nrow(predTable)
    }
  }
  # filter out rows with eGenes connected to un-expressed genes, write threshold table to output
  v = filter(v, !is.na(prediction.rate.inEnhancer))
  fName = file.path(outDir, "thresholdTables", method, paste0(tissue, ".", biosample, ".tsv"))
  write.table(v, file=fName, sep="\t", quote=F, row.names=F, col.names=T)

