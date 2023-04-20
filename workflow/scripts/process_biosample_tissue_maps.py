import pandas as pd
import os

#methodsFile = "/oak/stanford/groups/engreitz/Users/sheth/eQTLEnrichment-thresholding/eQTLEnrichment/config/config-baseline-predictors-newKeys.tsv"
#methods_config = pd.read_csv(methodsFile, sep='\t')

biosampleKeys = methods_config['sampleKey']
methods = methods_config['method']

biosample_maps = []
GTExTissue_maps = []
for i in range(len(methods)):
  if biosampleKeys[i]=="None":
    biosample_map.append("")
    GTExTissue_map.append("")
  else:
    biosample_map_this = []
    GTExTissue_map_this = []
    key = pd.read_csv(biosampleKeys[i], sep='\t')
    for i in range(len(key["biosample"])):
      if key['GTExTissue'][i]!="":
        biosample_map_this = biosample_map_this + key["biosample"][i]
        GTExTissue_map_this = GTExTissue_map_this + key["GTExTissue"][i]

    biosample_maps.append(biosample_map_this)
    GTExTissue_maps.append(GTExTissue_map_this)

methods_config['biosample_map'] = samples # add samples to config
methods_config['GTExTissue_map'] = files 

print(methods_config['biosample_map'])
print(methods_config['GTExTissue_map'])
