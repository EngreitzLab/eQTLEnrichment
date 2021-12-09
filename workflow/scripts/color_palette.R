suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(colorspace)})

main <- function() {
  ## get files from snakemake
  # methosd
  # colors
  user.inputs = (snakemake@params$user_inputs)
  print(user.inputs)
  methods = (snakemake@params$names) %>% strsplit(" ") %>% unlist
  print(methods)
  outFile = (snakemake@output$colorPalette)
  
  cp = data.frame(method=methods, hex="")
  count = 0;
  for (i in 1:nrow(cp)){
    if (user.inputs[i]!="None"){
      cp$hex[i]=user.inputs[i]
    } else {
    count=count+1
    }
  }
  

  # generate colorspace palette and fill in 
  cols = qualitative_hcl(n=count, palette="Set 2"); 
  count=1;
  for (i in 1:nrow(cp)){
    if (cp$hex[i]==""){
      cp$hex[i] = cols[count]
      count=count+1
    }
  }
  
  # add standards
  #standards = data.frame(method=c("TSS within 100 kb", "Closest TSS", "Closest gene body", "Proximity"),
  #                       hex = c("#7A7A7A", "#A8A8A8", "#CDCDCD", "#E2E2E2"))
  #cp = rbind(cp, standards)
  #cpList = split(f=cp$method, x=cp$hex)
  

	# write table
	saveRDS(cpList, file=outFile)
	
}

main()
