import pandas as pd
import os
import numpy as np

# for testing
#methodsFile = "/oak/stanford/groups/engreitz/Users/sheth/eQTLEnrichment-thresholding/eQTLEnrichment/config/config-baseline-predictors-newKeys.tsv"
#methods_config = pd.read_csv(methodsFile, sep='\t')

# read in parameters from config
nSteps = methods_config['nSteps']
maxScore = methods_config['maxScore']
minScore = methods_config['minScore']
inverse_predictor = methods_config['inverse_predictor']

thresholdSpans = [] # initialize list

for i in range(len(nSteps)): # iterate through number of methods
  maxScore_this = maxScore[i]
  minScore_this = minScore[i]
  nSteps_this = nSteps[i]
  inverse_this = inverse_predictor[i]

  if (maxScore_this-minScore_this)/(nSteps_this-1)>=nSteps_this: # if values of span will be integers
    span = np.linspace(minScore_this, maxScore_this, nSteps_this, dtype=int)
    span = span.tolist()

  else: # if values of span will be decimals
    span = np.linspace(minScore_this, maxScore_this, nSteps_this)
    span = [round(elem, ndigits=4) for elem in span] # converts to list

  if inverse_this==True:
    span = [i * -1 for i in span]

  thresholdSpans.append(span) # add that method to the list

methods_config['thresholdSpan'] = thresholdSpans # add to config
#print(methods_config['thresholdSpan'])
