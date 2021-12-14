import pandas as pd
import os

biosampleKeys = methods_config['sampleKey']
methods = methods_config['method']

samples = []
for i in range(len(methods)):
  if biosampleKeys[i]=="None":
    samples.append("")
  else:
    key = pd.read_csv(biosampleKeys[i], sep='\t')
    sampleList = key["biosample"].tolist()
    samples.append(sampleList)

methods_config['biosamples'] = samples # add samples to config
