import pandas as pd
import os
import numpy as np

#methodsFile = "/oak/stanford/groups/engreitz/Users/sheth/eQTLEnrichment-thresholding/eQTLEnrichment/config/config-baseline-predictors-newKeys.tsv"
#methods_config = pd.read_csv(methodsFile, sep='\t')

nSteps = methods_config['nSteps']
maxScore = methods_config['maxScore']

thresholdSpans = []

for i in range(len(nSteps)):
  maxScore.this = maxScore[i]
  nSteps.this = nSteps[i]
  if maxScore.this/(nSteps.this-1)>=nSteps.this:
    span = np.linspace(0, maxScore.this, nSteps.this, dtype=int)
    span = span.tolist()
  else:
    span = np.linspace(0, maxScore.this, nSteps.this)
    span = [round(elem, ndigits=4) for elem in span] # converts to list
    
  thresholdSpans.append(span)

methods_config['thresholdSpan'] = thresholdSpans # add to config
#print(methods_config['thresholdSpan'])
