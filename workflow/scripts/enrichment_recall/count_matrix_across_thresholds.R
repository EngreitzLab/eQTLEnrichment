suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(data.table)})

# iterate through each biosample intersection file indiviually and generate count matrix
## get files from snakemake
method = (snakemake@wildcards$method)
varPredInt_files = (snakemake@input$varPredInt) %>% strsplit(" ") %>% unlist() 
thresholdSpan = read.table(snakemake@input$thresholdSpan, sep="\t", header=FALSE) %>% setNames("threshold")
biosamples = (snakemake@params$biosamples) %>% strsplit(" ") %>% unlist()
varPerTissue = (snakemake@input$variantsPerTissue)
outFile = (snakemake@output$countMatrix)

# intialize df with cols = tissue, rows = biosample
varPerTissue = fread(varPerTissue, sep="\t", header=TRUE)
GTExTissues = varPerTissue$tissue

## read in and filter variant intersection
for (b in 1:length(biosamples)){
	sample.this = biosamples[b]
    varIntFile = varPredInt_files[b]
    
    size.file = file.info(varIntFile)$size
    size.threshold = 100 # in bytes, for empty file
	
	counts = data.frame(matrix(ncol=length(GTExTissues), nrow=length(thresholdSpan$threshold)))
	colnames(counts) = GTExTissues
	counts$Biosample = sample.this
	counts$threshold = thresholdSpan$threshold

	if (size.file<size.threshold){
     		 for (tissue in GTExTissues)
			 	counts[[tissue]] = 0
    } else { # read in file
		variantsInt = read.table(varIntFile, header=FALSE, sep="\t") %>% 
						setNames(c("varChr", "varStart", "varEnd", "variantID", "gene", "GTExTissue", "PIP", "distance_bin",
							"enhChr", "enhStart", "enhEnd", "Biosample", "TargetGene", "score"))
		for (t in 1:nrow(thresholdSpan)){ 	# iterate through thresholds
			score.thresh = thresholdSpan$threshold[t]
			variantsThresh = dplyr::filter(variantsInt, score>=score.thresh)  %>% dplyr::select(varChr, varStart, varEnd, GTExTissue, Biosample)
			for (tissue in GTExTissues){
				variantsTissue = dplyr::filter(variantsThresh, GTExTissue==tissue, Biosample==sample.this) %>% distinct()
				counts[counts$threshold==score.thresh,tissue] = nrow(variantsTissue)
      		}
    	}
	}
	if (b==1){
		counts_all = counts
	} else {
		counts_all = rbind(counts_all, counts)
  }
}

write.table(counts_all, file=outFile, sep="\t", quote=F, row.names=F, col.names=T)


