import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.stats import mannwhitneyu
from stats_util import is_outlier
from random import sample
import pandas as pd
import numpy as np
import sys

score = sys.argv[1] #'gwRVIS' #'JARVIS'
region = sys.argv[2] #'intron'



df = pd.read_csv('bed/JARVIS.prediction_scores.bed', sep='\t')
df.columns = ['chr', 'start', 'end', 'genomic_class', 'gwRVIS', 'JARVIS']
print(df.head())



jarvis_df = df[ ['genomic_class', score] ]
print(jarvis_df.head())

jarvis_df = jarvis_df.pivot(columns='genomic_class')
jarvis_df = jarvis_df[score].reset_index()
jarvis_df.columns.name = None
jarvis_df.drop(['index'], axis=1, inplace=True)
print(jarvis_df.head())


fig, ax = plt.subplots(figsize=(12, 12))
legend_dict = {'ccds': '#1f78b4',  'intergenic': '#e31a1c', 'intron': '#b2df8a', 'lincrna': '#737373', 'ucne': '#a6cee3', 'utr': '#33a02c', 'vista': '#6a3d9a', 'SV-intergenic': 'black', 'SV-intron': 'black', 'SV-utr': 'black' }


sv_jarvis_df = pd.read_csv('../../gwRVIS/other_datasets/structural-variants/jarvis-sv-bed/jarvis.SV.' + region + '.bed', sep='\t', header=None)
sv_jarvis_df.columns = ['AF', 'gwRVIS', 'JARVIS']
print(sv_jarvis_df.head())
print(sv_jarvis_df.shape)
sv_jarvis_df = sv_jarvis_df.loc[ sv_jarvis_df.AF >= 0.001, :]
print('Keeping common only (MAF > 0.1%):', sv_jarvis_df.shape)


bw=0.6
sv_scores = sv_jarvis_df[score].copy()
if score == 'gwRVIS':
	sv_scores = sv_scores[ ~is_outlier(sv_scores, 3.5) ]
sv_scores.plot.kde(bw_method=bw, color=legend_dict['SV-' + region])	

print('SV-' + region + ' median:', np.median(sv_scores))

if score == 'JARVIS':
	print('SV-' + region + ' mean:', np.mean(sv_scores))
	high_scores = sv_scores[ sv_scores > 0.6 ]
	low_scores = sv_scores[ sv_scores < 0.4 ]
	print('SV-' + region + ' mean (high regions > 0.6):', np.mean(high_scores))
	print('SV-' + region + ' mean (low regions < 0.4):', np.mean(low_scores))


res = mannwhitneyu(jarvis_df[region].sample(n=sv_jarvis_df.shape[0]), sv_jarvis_df[score])
print('>> [Mann-Whitney U] ' + score + ': ' + str(res.statistic) + ', P-value: ' + str(res.pvalue) + "\r\n")


if region == 'all':
	regions = jarvis_df.columns
else:
	regions = [region]


for col in regions: 
	print(">", col)
	col_data = jarvis_df[col].dropna()
	if score == 'gwRVIS':
		col_data = col_data[ ~is_outlier(col_data, 3.5) ]

	print('Whole region ' + region + ' median:', np.median(col_data))
	if score == 'JARVIS':
		print('Whole region ' + region + ' mean:', np.mean(col_data))
		high_scores = col_data[ col_data > 0.6 ]
		low_scores = col_data[ col_data < 0.4 ]
		print('Whole region ' + region + ' mean (high regions > 0.6):', np.mean(high_scores))
		print('Whole region ' + region + ' mean (low regions < 0.4):', np.mean(low_scores))

	col_data.plot.kde(bw_method=bw, color=legend_dict[col])	


if region != 'all':
	regions.append('SV-' + region)

patchList = []
for key in regions:
	data_key = mpatches.Patch(color=legend_dict[key], label=key)
	patchList.append(data_key)
plt.legend(handles=patchList)

fig.savefig('figs/' + score + '_distribution.' + region + '.pdf', bbox_inches='tight')
