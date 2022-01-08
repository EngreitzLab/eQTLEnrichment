import pandas as pd
import os
import numpy as np

#methodsFile = "/oak/stanford/groups/engreitz/Users/sheth/eQTLEnrichment-thresholding/eQTLEnrichment/config/config-baseline-predictors-newKeys.tsv"
#methods_config = pd.read_csv(methodsFile, sep='\t')

#nSteps = methods_config['nSteps']
#maxScore = methods_config['maxScore']

nSteps = 10
maxScore = 1000

span = linspace(0, maxScore, nSteps)
# if maxScore/nSteps is greater than 1...
# create vector and round each value to next greatest integer

# if maxScore/nSteps is less than 1...
# create vector and truncate each value to two sig-figs??

print(span)

#methods_config['thresholdSpan'] = span # add span to config

#print(methods_config['thresholdSpan])
