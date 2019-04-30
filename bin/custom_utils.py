import os, sys
from subprocess import call
import pandas as pd
from pathlib import Path
import yaml
import re


def create_out_dir(config_file):
	"""
	Create out directory for current run 
	based on run parameters - if it doesn't exist

	Input:
		config_file: (str)

	Output:
		out_dir: output directory path
	"""

	config_params = get_config_params(config_file)
	win_len = config_params['win_len']
	MAF_thres = config_params['MAF_thres']
	filter_outliers_before_regression = config_params['filter_outliers_before_regression']
	dataset = config_params['dataset']
	run_identifier = config_params['run_identifier']
	
	variants_table_dir = config_params['variants_table_dir']
	variant_and_popul_id = re.sub(".*filtered_variant_tables-", "", variants_table_dir)

	base_out_dir = '../out'
	out_dir = base_out_dir + '/' + dataset + '-' + run_identifier + '-winlen_' + str(win_len) + '.MAF_' + str(MAF_thres) + '.' + variant_and_popul_id

	if not os.path.exists(base_out_dir):
		os.makedirs(base_out_dir, exist_ok=True)
	if not os.path.exists(out_dir):     
		os.makedirs(out_dir, exist_ok=True)

	if not os.path.exists(out_dir + '/' + config_file):
		call(['cp', config_file, out_dir])

	return out_dir


def get_config_params(config_file):
	"""
	Read config file and store run parameters 
	in the 'config_params' dictionary

	Input:
		config_file: (str)

	Output: config_params: (dict)
	"""

	config_params = {}
	config_file = Path(config_file)
	
	with open(config_file, 'r') as yml_file:
		conf = yaml.load(yml_file)
		
	for group, _ in conf.items():
		for param, val in conf[group].items():
			config_params[param] = val
	
	return config_params


def set_intersection_checker(file1, file2):

	with open(file1) as f:
		content = f.readlines()
	set1 = [x.strip() for x in content]
	set1 = set(set1)

	with open(file2) as f:
		content = f.readlines()
	set2 = [x.strip() for x in content]
	set2 = set(set2)

	intersect = set.intersection(set1, set2)

	print(intersect)
	print('Number of intersection elements:', len(intersect))



# splitDataFrameList() adapted from:
# https://gist.github.com/jlln/338b4b0b55bd6984f883
def splitDataFrameList(df, target_column, separator):
    ''' 
    Input:
    	df = dataframe to split,
    	
	target_column = the column containing the values to split
    	
	separator = the symbol used to perform the split
    
    Output: 
	a dataframe with each entry for the target column separated, with each element moved into a new row. 
        The values in the other columns are duplicated across the newly divided rows.
    '''
    def splitListToRows(row, row_accumulator, target_column, separator):
        split_row = row[target_column].split(separator)
        for s in split_row:
            new_row = row.to_dict()
            new_row[target_column] = s
            row_accumulator.append(new_row)
    new_rows = []
    df.apply(splitListToRows,axis=1,args = (new_rows, target_column, separator))
    new_df = pd.DataFrame(new_rows)
    return new_df


if __name__ == '__main__':

	config_file = sys.argv[1]
	config_params = get_config_params(config_file)
	print(config_params)
	
