import os
import pandas as pd

MAX_MEM_MB = 250 * 1000  # 250GB

def determine_mem_mb(wildcards, input, attempt, min_gb=8):
	# Memory resource calculator for snakemake rules
	input_size_mb = input.size_mb
	if ".gz" in str(input):
		input_size_mb *= 8  # assume gz compressesed the file <= 8x
	attempt_multiplier = 2 ** (attempt - 1)  # Double memory for each retry
	mem_to_use_mb = attempt_multiplier *  max(4 * input_size_mb, min_gb * 1000)
	return min(mem_to_use_mb, MAX_MEM_MB)

def add_biosamples_and_files_to_config():
	biosampleKeys = methods_config['sampleKey']
	methods = methods_config['method']

	samples = []
	files = []
	for i in range(len(methods)):
		if biosampleKeys.iloc[i]=="None":
			samples.append("")
			files.append("")
		else:
			key = pd.read_csv(biosampleKeys.iloc[i], sep='\t')
			sampleList = key["biosample"].tolist()
			fileList =  key["predictionFile"].tolist()
			fileDict = dict(zip(sampleList, fileList))
			samples.append(sampleList)
			files.append(fileDict)

	methods_config['biosamples'] = samples # add samples to config
	methods_config['predFiles'] = files 

def process_biosample_tissue_maps():
	biosampleKeys = methods_config['sampleKey']
	methods = methods_config['method']

	biosample_maps = []
	GTExTissue_maps = []
	for i in range(len(methods)):
		tissues_this = []
		biosamples_this = []
		if not biosampleKeys.iloc[i]=="None":
			key = pd.read_csv(biosampleKeys.iloc[i], sep='\t')
			key = key.dropna(subset='GTExTissue')
			# make lists, accounting for replicate matches
			for index, row in key.iterrows():
				tissues_this_row = [x.strip() for x in row["GTExTissue"].split(',')]
				biosamples_this_row = [row['biosample'].strip() for x in tissues_this_row]
				tissues_this.extend(tissues_this_row)
				biosamples_this.extend(biosamples_this_row)

		biosample_maps.append(biosamples_this)
		GTExTissue_maps.append(tissues_this)

	methods_config['GTExTissue_map'] = GTExTissue_maps
	methods_config['biosample_map'] = biosample_maps


def flatten(list_of_lists):
	return [x for xs in list_of_lists for x in xs]