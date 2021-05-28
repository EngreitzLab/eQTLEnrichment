suppressPackageStartupMessages({
  library(biomaRt)
  library(dplyr)
  library(tidyr)
  library(optparse)
  library(stringr)
  library(gdata)})


main <- function() {
  this.biosample = (snakemake@wildcards$biosample)
  sampleKeyFile = (snakemake@input$sampleKey)
  method.geneUniverse.file = (snakemake@input$geneUniverse)
  method.TSS.file = (snakemake@params$TSS)
  GTEx.geneUniverse.file = (snakemake@params$GTExGeneUniverse)

  ## decide on TSS file to use
   if (is.na(sampleKeyFile)) {
     TSS.file =  method.TSS.file
     } else {
        sampleKey = read.table(sampleKeyFile, header=TRUE, fill=TRUE, row.names=NULL)
        if ("TSSFile" %in% colnames(sampleKey)){
        filteredKey = filter(sampleKey, biosample==this.biosample)
        this.TSS = filteredKey$TSSFile[1]
        # if blank, use default; if not
        if (is.na(this.TSS) || this.TSS=="None") {
          TSS.file =  method.TSS.file
        } else {
          TSS.file = this.TSS
        }
      } else {
         TSS.file = method.TSS.file
      }
    }
 
  ## decide on gene file to use
   if (is.na(sampleKeyFile)) {
     gene.file = method.geneUniverse.file
   } else {
     sampleKey = read.table(sampleKeyFile, header=TRUE, sep="\t", fill=TRUE)
     if ("geneFile" %in% colnames(sampleKey)){
       filteredKey = filter(sampleKey, biosample==this.biosample)
       this.genes = filteredKey$geneFile[1]
       # if blank, use default; if not
       if (is.na(this.genes) || this.genes=="None") {
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
  
  # filter
  TSS = filter(TSS, gene %in% GTEx.geneUniverse$gene)
  genes = filter(genes, gene %in% GTEx.geneUniverse$gene)

  write.table(TSS, file=snakemake@output$specificTSSFile, sep="\t", quote=F, row.names=F, col.names=F)
  write.table(genes, file=snakemake@output$specificGeneFile, sep="\t", quote=F, row.names=F, col.names=F)
  
}

main()
