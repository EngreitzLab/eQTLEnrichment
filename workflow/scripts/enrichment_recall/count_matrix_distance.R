suppressPackageStartupMessages({library(dplyr)
  library(tidyr)
  library(stringr)
  library(data.table)})

# iterate through each biosample intersection file indiviually and generate count matrix
# "clear" from memory at each iteration

## get files from snakemake
score.thresh = (snakemake@params$threshold) %>% as.numeric()
distance.this = (snakemake@wildcards$distance_max) %>% as.numeric()
biosamples = (snakemake@params$biosamples) %>% strsplit(" ") %>% unlist()
varPredInt_files = (snakemake@input$varPredInt) %>% strsplit(" ") %>% unlist() 
varPerTissue = (snakemake@input$variantsPerTissue)
outFile = (snakemake@output$countMatrix)

# initialize df wtih columns = tissues, rows = biosamples
varPerTissue = fread(varPerTissue, sep="\t", header=TRUE)
GTExTissues = varPerTissue$tissue

counts = data.frame(matrix(ncol=length(GTExTissues), nrow=length(biosamples)))
colnames(counts) = GTExTissues
counts$Biosample = biosamples

## read in and filter variant intersection
# columns: 1-3 (loc), 4 (hgID), 5 (tissue), 6 (ens_id), 7 (PIP), 8 (TPM), 9 (distance), 
# 10-12 (enhancer loc), 13 (enhancer cell type), 14 (enhancer target gene hgnc), 15 (enhancer score)
for (i in 1:length(biosamples)){
  sample.this = biosamples[i]
  varIntFile = varPredInt_files[i]
  
  size.file = file.info(varIntFile)$size
  size.threshold = 100 # in bytes, for empty file
  
  # if intersection file is empty, set counts for this biosample across all tissues to 0
  if (size.file<size.threshold){
    for (tissue in GTExTissues)
      counts[counts$Biosample==sample.this,tissue] = 0   
    
  } else { # otherwise, read in file
    variantsInt = read.table(varIntFile, header=FALSE, sep="\t") %>% 
      setNames(c("varChr", "varStart", "varEnd", "variantID", "gene", "GTExTissue", "PIP", "distance_bin",
                 "enhChr", "enhStart", "enhEnd", "Biosample", "TargetGene", "score"))
    # filter to distance threshold and select columns
    if (distance.this==30000){ # "all variants" category
      variantsInt = dplyr::select(variantsInt, variantID, GTExTissue, Biosample)
    } else{
      variantsInt = dplyr::filter(variantsInt, score>=score.thresh, distance_bin==distance.this) %>% 
        dplyr::select(variantID, GTExTissue, Biosample)
    }

    for (tissue in GTExTissues){
      variantsTissue = dplyr::filter(variantsInt, GTExTissue==tissue, Biosample==sample.this) %>% distinct()
     counts[counts$Biosample==sample.this,tissue] = nrow(variantsTissue)
    }
    rm(variantsInt)
  }
}

write.table(counts, file=outFile, sep="\t", quote=F, row.names=F, col.names=T)
    

