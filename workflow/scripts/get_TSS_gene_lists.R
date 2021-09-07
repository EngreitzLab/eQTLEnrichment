suppressPackageStartupMessages({
  library(biomaRt)
  library(dplyr)
  library(tidyr)
  library(optparse)
  library(stringr)
  library(gdata)})


main <- function() {
  this.biosample = (snakemake@wildcards$biosample)
  sampleKeyFile = (snakemake@params$sampleKey)
  method.geneUniverse.file = (snakemake@input$geneUniverse)
  method.TSS.file = (snakemake@params$TSS)
  GTEx.geneUniverse.file = (snakemake@params$GTExGeneUniverse)

  ## decide on TSS file to use
   if (is.na(sampleKeyFile) || sampleKeyFile=="None") {
     TSS.file =  method.TSS.file
     } else {
        sampleKey = read.table(sampleKeyFile, header=TRUE, fill=TRUE, row.names=NULL, sep='\t')
        if ("TSSFile" %in% colnames(sampleKey)){
        filteredKey = filter(sampleKey, biosample==this.biosample)
        this.TSS = filteredKey$TSSFile[1]
        # if blank, use default; if not
        if (is.na(this.TSS) || this.TSS=="None" || this.TSS=="") {
          TSS.file =  method.TSS.file
        } else {
          TSS.file = this.TSS
        }
      } else {
         TSS.file = method.TSS.file
      }
    }
 
  ## decide on gene file to use
   if (is.na(sampleKeyFile) || sampleKeyFile=="None") {
     gene.file = method.geneUniverse.file
   } else {
     sampleKey = read.table(sampleKeyFile, header=TRUE, row.names=NULL, fill=TRUE, sep='\t')
     if ("geneFile" %in% colnames(sampleKey)){
       filteredKey = filter(sampleKey, biosample==this.biosample)
       this.genes = filteredKey$geneFile[1]
       # if blank, use default; if not
       if (is.na(this.genes) || this.genes=="None" || this.genes=="") {
          gene.file = method.geneUniverse.file
       } else {
         gene.file = this.genes
       }
     } else {
       gene.file = method.geneUniverse.file
     }
   }     

  ## filter genes and TSS (specific to method/biosample) to gene universe for this method (by GTEx gene universe)
  # read in files
  GTEx.geneUniverse = read.table(file=GTEx.geneUniverse.file, header=FALSE)
  colnames(GTEx.geneUniverse) = c("chr", "start", "end", "gene", "score", "strand")
  TSS = read.table(file=TSS.file, header=FALSE)
  colnames(TSS) = c("chr", "start", "end", "gene", "score", "strand")
  genes = read.table(file=gene.file, header=FALSE)
  colnames(genes) = c("chr", "start", "end", "gene", "score", "strand")
  
  # make sure gene column of TSS is actually gene name
  for (i in 1:nrow(TSS)){
    y = str_split(TSS$gene[i], "-") %>% unlist()
    if (length(y)>1){
      z = str_split(y[2], "_") %>% unlist()
      TSS$gene[i] = z[2]
    }
  }

  # make sure gene column of gene list is actually gene name
  for (i in 1:nrow(genes)){
    y = str_split(genes$gene[i], "-") %>% unlist()
    if (length(y)>1){
      genes$gene[i] = y[1]
    }
  }

  # filter
  #TSS = dplyr::filter(TSS, is.element(gene,GTEx.geneUniverse$gene))
  TSS = TSS[is.element(TSS$gene, GTEx.geneUniverse$gene), ]
  #genes = dplyr::filter(genes, iselement(gene, GTEx.geneUniverse$gene))
  genes = genes[is.element(genes$gene, GTEx.geneUniverse$gene), ]
  
  write.table(TSS, file=snakemake@output$specificTSSFile, sep="\t", quote=F, row.names=F, col.names=F)
  write.table(genes, file=snakemake@output$specificGeneFile, sep="\t", quote=F, row.names=F, col.names=F)
  
}

main()
