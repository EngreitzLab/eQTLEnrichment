import pandas as pd
import os

#methodsFile = "/oak/stanford/groups/engreitz/Users/sheth/eQTLEnrichment-byBiosample/eQTLEnrichment/config/config-baseline-predictors-newKeys.tsv"
#methods_config = pd.read_csv(methodsFile, sep='\t')

biosampleKeys = methods_config['sampleKey']
methods = methods_config['method']

samples = []
files = []
for i in range(len(methods)):
  if biosampleKeys[i]=="None":
    samples.append("")
    files.append("")
  else:
    key = pd.read_csv(biosampleKeys[i], sep='\t')
    sampleList = key["biosample"].tolist()at
    fileList =  key["predictionFile"].tolist()
    fileDict = dict(zip(sampleList, fileList))
    samples.append(sampleList)
    files.append(fileDict)

methods_config['biosamples'] = samples # add samples to config
methods_config['predFiles'] = files 

#print(methods_config['predFiles'])
