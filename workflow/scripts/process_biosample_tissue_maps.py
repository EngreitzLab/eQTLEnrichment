import pandas as pd
import os
import math

biosampleKeys = methods_config['sampleKey']
methods = methods_config['method']

biosample_maps = []
GTExTissue_maps = []
for i in range(len(methods)):
  if biosampleKeys[i]=="None":
    biosample_map_this = []
    GTExTissue_map = []
  else:
    key = pd.read_csv(biosampleKeys[i], sep='\t')

    key = key.dropna()
    biosample_map_this = key["biosample"].tolist()
    GTExTissue_map_this = key["GTExTissue"].tolist()

  biosample_maps.append(biosample_map_this)
  GTExTissue_maps.append(GTExTissue_map_this)


methods_config['biosample_map'] = samples # add samples to config
methods_config['GTExTissue_map'] = files 
