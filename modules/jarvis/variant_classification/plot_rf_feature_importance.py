import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys


def plot_feat_imp_per_class(df, genomic_class):
	fig, ax = plt.subplots(figsize=(8, 12))

	ax.barh(df['feature'], df['imp_score'], height=0.6)
	plt.title(genomic_class)
	fig.savefig(genomic_class + "_feat_imp.pdf", bbox_inches='tight')



def merge_dfs(all_dfs):

	rand_key = [*all_dfs][0]
	merged_df = all_dfs[rand_key]

	del all_dfs[rand_key]

	print(merged_df.head())
	print(all_dfs.keys())

	for key in all_dfs.keys():
		merged_df = merged_df.merge(all_dfs[key], left_on='feature', right_on='feature', how='outer')

		
	merged_df.index = merged_df['feature']
	del merged_df['feature']

	# min-max normalisation
	#merged_df = (merged_df - merged_df.min()) / (merged_df.max() - merged_df.min())


	merged_df['avg'] = merged_df.mean(axis=1)

	merged_df.reset_index(inplace=True)
	merged_df.rename(columns={'avg': 'imp_score'}, inplace=True)
	merged_df.sort_values(by='imp_score', inplace=True)
	print(merged_df.head())


	
	
	plot_feat_imp_per_class(merged_df, 'avg_feat_imp')



if __name__ == "__main__":

	genomic_classes = ['intergenic', 'lincRNA', 'utr', 'intergenic-utr-lincrna-ucne-vista']
	all_dfs = {}


	for genomic_class in genomic_classes:
		print(">" + genomic_class)

		df = pd.read_csv(genomic_class + "_RF_feat_imp.with_dupl.txt", sep='\t', header=None)
		df.columns = ['feature', 'imp_score']
		plot_feat_imp_per_class(df, genomic_class)

		all_dfs[genomic_class] = df


	# get average feature importance from all cases
	merge_dfs(all_dfs)
