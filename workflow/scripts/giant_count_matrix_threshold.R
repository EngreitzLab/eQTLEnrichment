suppressPackageStartupMessages({library(dplyr)
  library(tidyr)
  library(stringr)
  library(data.table)})

# iterate through each biosample intersection file indiviually and generate count matrix
## get files from snakemake
method = (snakemake@wildcards$method)
#threshold = (snakemake@params$threshold)
thresholdSpan = read.table(snakemake@input$thresholdSpan, sep="\t", header=FALSE) %>% setNames("threshold")
print(thresholdSpan)
biosamples = (snakemake@params$biosamples) %>% strsplit(" ") %>% unlist() %>% sort()
GTExTissues = (snakemake@params$GTExTissues) %>% strsplit(" ") %>% unlist() %>% sort()
outDir = (snakemake@params$outDir)
outFile = (snakemake@output$countMatrix)

counts = data.frame(matrix(ncol=length(GTExTissues), nrow=length(biosamples)))
colnames(counts) = GTExTissues
counts$Biosample = biosamples

## read in and filter variant intersection
# columns: 1-3 (loc), 4 (hgID), 5 (tissue), 6 (ens_id), 7 (PIP), 8 (TPM), 9 (distance), 
# 10-12 (enhancer loc), 13 (enhancer cell type), 14 (enhancer target gene hgnc), 15 (enhancer score)

for (j in 1:nrow(thresholdSpan)){
  score.thresh=thresholdSpan$threshold[j]
  print(score.thresh)
  
  for (i in 1:length(biosamples)){
    sample.this = biosamples[i]
    varIntFile = file.path(outDir, method, sample.this, "GTExVariants-enhancerPredictionsInt.tsv.gz")
    
    size.file = file.info(varIntFile)$size
    size.threshold = 100 # in bytes, for empty file
    
    # if intersection file is empty, set counts for this biosample across all tissues to 0
    if (size.file<size.threshold){
      for (tissue in GTExTissues)
        counts[counts$Biosample==sample.this, tissue] = 0   
      
    } else { # otherwise, read in file
      variantsInt = read.table(varIntFile, header=FALSE, sep="\t") %>% 
        setNames(c("varChr", "varStart", "varEnd", "hgID", "GTExTissue", "gene", "PIP", "TPM", "distance",
                   "enhChr", "enhStart", "enhEnd", "Biosample", "TargetGene", "score"))
      
      # filter to score threshold and select columns
      variantsInt = dplyr::filter(variantsInt, score>=score.thresh) %>% 
        dplyr::select(varChr, varStart, varEnd, GTExTissue, Biosample)
      for (tissue in GTExTissues){
        variantsTissue = dplyr::filter(variantsInt, GTExTissue==tissue, Biosample==sample.this) %>% distinct()
        counts[counts$Biosample==sample.this,tissue] = nrow(variantsTissue)
      }
      counts$threshold = score.thresh
    }
  }
  if (j==1){
    counts_all = counts
  } else {
    counts_all = rbind(counts_all, counts)
  }
}

write.table(counts_all, file=outFile, sep="\t", quote=F, row.names=F, col.names=T)


