suppressPackageStartupMessages({library(dplyr)
  library(tidyr)
  library(stringr)})

# iterate through each biosample intersection file indiviually and generate count matrix
# "clear" from memory at each iteration

  main <- function() {
  ## get files from snakemake
  #dfFile = (snakemake@input$variantsPredictionsInt)
  method = (snakemake@wildcards$method)
  biosampleFile = (snakemake@input$samples)
  GTExTissues = (snakemake@params$GTExTissues) %>% strsplit(" ") %>% unlist() %>% sort()
  outDir = (snakemake@params$outDir)
  outFile = (snakemake@output$countMatrix)
  

	## read in files
  biosamples = read.table(biosampleFile, header=FALSE, stringsAsFactors=FALSE); biosamples = biosamples[[1]]
  counts = data.frame(matrix(ncol=length(GTExTissues), nrow=length(biosamples)))
  colnames(counts)=GTExTissues
  counts$Biosample=biosamples
  
  for (i in 1:length(biosamples)){
    sample.this = biosamples[i]
    int.file = file.path(outDir, method, sample.this, "GTExVariants-enhancerPredictionsInt.tsv.gz")
    size.file = file.info(int.file)$size
    size.threshold = 100 # in bytes, for empty file
    
    if (size.file<size.threshold){
      for (tissue in GTExTissues)
        counts[counts$Biosample==sample.this,tissue] = 0   
    } else {
    df = read.table(gzfile(int.file), header=FALSE)
    colnames(df) = c('var.chr','var.start','var.end','hgid','GTExTissue','eGene','PIP', 'TPM', 
                     'enh.chr', 'enh.start', 'enh.end', 'Biosample', 'TargetGene', 'Score')
    df = dplyr::select(df, hgid, GTExTissue, Biosample)
    
    # initialize counts matrix (ncol by 1 matrix for individual biosample)
    #counts.this = data.frame(matrix(ncol=length(GTExTissues), nrow=1))
    #counts.this$Biosample=sample.this
    
    # get counts
    for (tissue in GTExTissues){
      vars.tissue=filter(df, GTExTissue==tissue)
      vars.int = filter(vars.tissue, Biosample==sample.this) %>% distinct()
      counts[counts$Biosample==sample.this,tissue] = nrow(vars.int)
    }
    rm(int.file)
    }
  }
  
	# write matrix
	write.table(counts, file=outFile, sep="\t", quote=F, row.names=F, col.names=T)
}

main()
