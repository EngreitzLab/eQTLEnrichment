suppressPackageStartupMessages({library(dplyr)
  library(tidyr)
  library(stringr)})

# iterate through each biosample intersection file indiviually and generate count matrix
# "clear" from memory at each iteration

## get files from snakemake
method = (snakemake@wildcards$method)
score.thresh = (snakemake@input$threshold) %>% as.numeric()
distance.thresh = (snakemake@wildcards$distance) %>% as.numeric()
biosamples = (snakemake@params$biosamples) %>% strsplit(" ") %>% unlist() %>% sort()
GTExTissues = (snakemake@params$GTExTissues) %>% strsplit(" ") %>% unlist() %>% sort()
outDir = (snakemake@params$outDir)
outFile = (snakemake@output$countMatrix)

counts = data.frame(matrix(ncol=length(GTExTissues), nrow=length(biosamples)))
counts$Biosample = biosamples

## read in and filter variant intersection
# columns: 1-3 (loc), 4 (hgID), 5 (tissue), 6 (ens_id), 7 (PIP), 8 (TPM), 9 (distance), 
# 10-12 (enhancer loc), 13 (enhancer cell type), 14 (enhancer target gene hgnc), 15 (enhancer score)
for (i in 1:length(biosamples)){
  sample.this = biosamples[i]
  varIntFile = file.path(outDir, method, sample.this, "GTExVariants-enhancerPredictionsInt.tsv.gz")
  
  size.file = file.info(int.file)$size
  size.threshold = 100 # in bytes, for empty file
  
  # if intersection file is empty, set counts for this biosample across all tissues to 0
  if (size.file<size.threshold){
    for (tissue in GTExTissues)
      counts[counts$Biosample==sample.this,tissue] = 0   
    
  } else { # otherwise, read in file
    variantsInt = read.table(variantsIntFile, header=FALSE, sep="\t") %>% 
      setNames(c("varChr", "varStart", "varEnd", "hgID", "GTExTissue", "gene", "PIP", "TPM", "distance",
                 "enhChr", "enhStart", "enhEnd", "Biosample", "TargetGene", "score"))
    # filter to distance threshold and select columns
    variantsInt = dplyr::filter(variantsInt, score>=score.thresh, distance<=distance.thresh) %>% 
      dplyr::select(hgID, GTExTissue, biosample)
    for (tissue in GTExTissues){
      variantsTissue = dplyr::filter(variantsInt, GTExTissue==tissue, Biosample==sample.this) %>% distinct()
      counts[counts$Biosample==sample.this,tissue] = nrow(variantsTissue)
    }
    rm(variantsInt)
  }
}

write.table(counts, file=outFile, sep="\t", quote=F, row.names=F, col.names=T)
    

