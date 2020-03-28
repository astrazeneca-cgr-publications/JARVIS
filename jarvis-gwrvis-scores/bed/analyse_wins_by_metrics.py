import pandas as pd
import sys, os


def read_full_df():

	df = pd.read_csv('jarvis-phastcons-output/mean_scores_by_window.bed', sep='\t', header=None)
	df.columns = ['chr', 'start', 'end', 'gwrvis', 'jarvis', 'phastcons']
	print(df.head())


	df.sort_values(by='phastcons', ascending=True, inplace=True)
	print(df.head())
	print(df.shape)

	return df



if __name__ == '__main__':


	top_ratio = 0.01
	ccds_bed = 'ccds_regions_with_gene_names.merged.bed'

	full_df = read_full_df()

	least_conserv_df = full_df.head(int(top_ratio * full_df.shape[0]))
	print(least_conserv_df.shape) 


	least_conserv_most_intol_df = least_conserv_df.sort_values(by='gwrvis', ascending=True)


	least_conserv_most_intol_df.to_csv('./least_conserv_most_intol.csv', sep='\t', index=False, header=False)
