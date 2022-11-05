suppressPackageStartupMessages({library(dplyr)
  library(tidyr)
  library(stringr)})

# iterate through each biosample intersection file indiviually and generate count matrix
# "clear" from memory at each iteration

  main <- function() {
  ## get files from snakemake
  method = (snakemake@wildcards$method)
  threshold = (snakemake@wildcards$threshold)
  biosampleFile = (snakemake@input$samples)
  GTExTissues = (snakemake@params$GTExTissues) %>% strsplit(" ") %>% unlist() %>% sort()
  outDir = (snakemake@params$outDir)
  outFile = (snakemake@output$countMatrix)
  
	## read in files
  biosamples = read.table(biosampleFile, header=FALSE, stringsAsFactors=FALSE); biosamples = biosamples[[1]]
  counts = data.frame(matrix(ncol=length(GTExTissues), nrow=length(biosamples)))
  colnames(counts)=GTExTissues
  counts$Biosample=biosamples
  
  ## iterate through biosamples
  for (i in 1:length(biosamples)){
    sample.this = biosamples[i]
    int.file = file.path(outDir, method, sample.this, 
                         paste0("GTExVariants-enhancerPredictionsInt.", threshold, ".tsv.gz"))
    size.file = file.info(int.file)$size
    size.threshold = 100 # in bytes, for empty file
    
    # if intersection file is empty, set counts for this biosample across all tissues to 0
    if (size.file<size.threshold){
      for (tissue in GTExTissues)
        counts[counts$Biosample==sample.this,tissue] = 0   
      
    } else { # otherwise, read in file
    df = read.table(gzfile(int.file), header=FALSE)
    colnames(df) = c('var.chr','var.start','var.end','hgid','GTExTissue','eGene','PIP', 'TPM', 
                     'enh.chr', 'enh.start', 'enh.end', 'Biosample', 'TargetGene', 'Score')
    df = dplyr::select(df, hgid, GTExTissue, Biosample)
    
    # get counts for each GTEx tissue
    for (tissue in GTExTissues){
      vars.tissue=filter(df, GTExTissue==tissue)
      vars.int = filter(vars.tissue, Biosample==sample.this) %>% distinct()
      counts[counts$Biosample==sample.this,tissue] = nrow(vars.int)
    }
    
    # clear intersection file from memory
    rm(int.file)
    }
  }
  
	## write matrix to output file
	write.table(counts, file=outFile, sep="\t", quote=F, row.names=F, col.names=T)
}

main()
