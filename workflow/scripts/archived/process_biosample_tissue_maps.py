import pandas as pd
import os
import math

biosampleKeys = methods_config['sampleKey']
methods = methods_config['method']

biosample_maps = []
GTExTissue_maps = []
for i in range(len(methods)):
	if biosampleKeys[i]=="None":
		tissues_this = []
		biosamples_this = []
	else:
		key = pd.read_csv(biosampleKeys[i], sep='\t')
		key = key.dropna(subset='GTExTissue')
		# make lists, accounting for replicate matches
		for index, row in key.iterrows():
			tissues_this = [x.strip() for x in row["GTExTissue"].split(',')]
			biosamples_this = [row['biosample'].strip() for x in tissues_this]
		
	GTExTissue_maps.append(tissues_this)
	biosample_maps.append(biosamples_this)

# add to methods_config
methods_config['GTExTissue_map'] = GTExTissue_maps 
methods_config['biosample_map'] = biosample_maps 

