import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

import pandas as pd
import sys
import numpy as np
from scipy.stats import mannwhitneyu


def is_outlier(points, thresh=3.5):

	if len(points.shape) == 1:         
		points = points[:,None]

	median = np.median(points, axis=0)     

	diff = np.sum((points - median)**2, axis=-1)     

	diff = np.sqrt(diff)

	med_abs_deviation = np.median(diff)
	if len(points) == 1 and med_abs_deviation == 0.0:
		return np.array([False])

	modified_z_score = 0.6745 * diff / med_abs_deviation

	return modified_z_score > thresh

	


def plot_for_score_and_class(score, genomic_class):
	
	class_file = score + '.' + genomic_class + '.txt'
	if genomic_class == 'intron':
		sv_class_file = score + '.sv_intronic.txt'
	else:
		sv_class_file = score + '.sv_' + genomic_class + '.txt'
	

	# Gen-Pop
	print('\n', class_file)
	score_df = pd.read_csv(class_file, low_memory=False, header=None)
	print(score_df.shape[0])
	score_list = score_df.iloc[:, 0]
	del score_df
	print('Median:', np.median(score_list))
	print('Mean:', np.mean(score_list))

	# SV
	print(sv_class_file)
	sv_score_df = pd.read_csv(sv_class_file, low_memory=False, header=None)
	print(sv_score_df.shape[0])
	sv_score_list = sv_score_df.iloc[:, 0]
	del sv_score_df
	print('Median:', np.median(sv_score_list))
	print('Mean:', np.mean(sv_score_list))


	"""
	if score == 'gwrvis':
		print('Filtering outliers for gwrvis')
		score_list = score_list[ ~is_outlier(score_list)]
		sv_score_list = sv_score_list[ ~is_outlier(sv_score_list)]
		

	fig, ax = plt.subplots(figsize=(12, 9))

	bw_dict = {'utr': {'gwrvis': 0.5, 'jarvis': 0.1}, 
		   'intron': {'gwrvis': 0.5, 'jarvis': 0.1} }

	sns.kdeplot(score_list, bw=bw_dict[genomic_class][score], label=genomic_class, color='#3182bd')
	sns.kdeplot(sv_score_list, bw=bw_dict[genomic_class][score],  label='sv_' + genomic_class, color='#de2d26')

	fig.savefig(score + '.' + genomic_class + '.gen_pop_vs_sv_distr.pdf', bbox_inches='tight')
	"""

	res = mannwhitneyu(score_list, sv_score_list)
	print('Mann-Whitney U:', res)



if __name__ == '__main__':

	score = sys.argv[1]
	genomic_class = sys.argv[2]

	plot_for_score_and_class(score, genomic_class)
