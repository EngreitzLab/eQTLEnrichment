import pandas as pd
import os
import numpy as np

# for testing
#methodsFile = "/oak/stanford/groups/engreitz/Users/sheth/eQTLEnrichment-thresholding/eQTLEnrichment/config/config-baseline-predictors-newKeys.tsv"
#methods_config = pd.read_csv(methodsFile, sep='\t')

# read in parameters from config
nSteps = methods_config['nSteps']
maxScore = methods_config['maxScore']

thresholdSpans = [] # initialize list

for i in range(len(nSteps)): # iterate through number of methods
  maxScore.this = maxScore[i]
  nSteps.this = nSteps[i]
  if maxScore.this/(nSteps.this-1)>=nSteps.this: # if values of span will be integers
    span = np.linspace(0, maxScore.this, nSteps.this, dtype=int)
    span = span.tolist()
  else: # if values of span will be decimals
    span = np.linspace(0, maxScore.this, nSteps.this)
    span = [round(elem, ndigits=4) for elem in span] # converts to list
    
  thresholdSpans.append(span) # add that method to the list

methods_config['thresholdSpan'] = thresholdSpans # add to config
#print(methods_config['thresholdSpan'])
