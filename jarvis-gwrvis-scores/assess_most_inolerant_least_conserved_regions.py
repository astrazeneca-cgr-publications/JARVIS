import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
import sys


# Read phastCons primate (and then 46-way)
# Ref. file with phastCons-primate conservation scores (indexed with tabix):
# ../../other_datasets/conservation/phastCons46way_primates/bed/All_chr.phastCons46way.primates.sorted.bed.gz


# Intersect with gwRVIS


# Then select those that are least conserved (get average conervation score over each window?)


# Run Pathway analysis/GO enrichment


def plot_score_distribution_across_classes(df, score='gwRVIS', top_ratio=1, tail='intolerant'):

	fig, ax = plt.subplots(figsize=(12, 12))

	legend_dict = {'ccds': '#1f78b4',  'intergenic': '#e31a1c',
		       'intron': '#b2df8a', 'lincrna': '#737373',
		       'ucne': '#a6cee3', 'utr': '#33a02c', 'vista': '#6a3d9a' }

	print("\n(Temp) Plotting df:")
	print(df.head())
	print(df.tail())
	print(df.shape)

	for col in df.columns:
		print(">", col)

		col_data = df[col].dropna()
		if len(col_data) < 2:
			continue
		col_data.plot.kde(bw_method=1, color=legend_dict[col])	
		#col_data.plot.hist(color=legend_dict[col])	
		

	patchList = []
	for key in legend_dict:
		data_key = mpatches.Patch(color=legend_dict[key], label=key)
		patchList.append(data_key)
	plt.legend(handles=patchList)

	fig.savefig('figs/' + score + '_distribution.top_' + str(top_ratio) + 'perc.most_' + tail + '.pdf', bbox_inches='tight')



def plot_regions_by_intolerance(df, score='gwRVIS', top_ratio=1, tail='intolerant'):

	"""
	    Input df: has gwRVIS sorted in ascending order
	"""

	print("\n> Top-" + str(top_ratio) + "% most " + tail + ":")
	if tail == 'intolerant':
		sub_df = df.head( int(df.shape[0] * top_ratio / 100.0) )
	elif tail == 'tolerant':
		sub_df = df.tail( int(df.shape[0] * top_ratio / 100.0) )
	print(sub_df.shape)	
	
	# save sub-df into file
	sub_df.to_csv('top_' + str(top_ratio) + 'perc.most_' + tail + '.' + score + '.tsv', sep='\t', header=False, index=False)

	sub_df = sub_df[ ['genomic_class', score] ]

	sub_df = sub_df.pivot(columns='genomic_class')
	sub_df = sub_df[score].reset_index()
	sub_df.columns.name = None
	sub_df.drop(['index'], axis=1, inplace=True)

	plot_score_distribution_across_classes(sub_df, score=score, top_ratio=top_ratio, tail=tail)



if __name__ == '__main__':

	df = pd.read_csv('bed/JARVIS.prediction_scores.bed', sep='\t')
	df.columns = ['chr', 'start', 'end', 'genomic_class', 'gwRVIS', 'JARVIS']
	print(df.head(20))

	# Sort gwRVIS to focus on the top 1,5,10% most intolerant
	# chr   start     end genomic_class    gwRVIS        JARVIS
	gwrvis_df = df[ ['chr', 'start', 'end', 'genomic_class', 'gwRVIS'] ].copy()
	gwrvis_df.sort_values(by='gwRVIS', inplace=True, ascending=True)
	print(gwrvis_df.shape)
	print(gwrvis_df.head())

	plot_regions_by_intolerance(gwrvis_df, score='gwRVIS', top_ratio=0.1, tail='intolerant')
	plot_regions_by_intolerance(gwrvis_df, score='gwRVIS', top_ratio=0.1, tail='tolerant')

