import pandas as pd
import sys


def group_wins_in_chrom(chrom):

	df = pd.read_csv(base_dir + '/JARVIS_phastcons.full_table.chr' + str(chrom) + '.bed.gz', sep='\t', compression='gzip', header=None)

	df.columns = ['chrom', 'start', 'end', 'gwRVIS', 'JARVIS', 'chr_phastC', 'start_phastC', 'end_phastC', 'phastCons', 'overlap']
	df.drop(['chr_phastC', 'start_phastC', 'end_phastC', 'overlap'], axis=1, inplace=True)
	print(df.head())

	# group by window coordinates and get mean scores
	group_df = df.groupby(['chrom', 'start', 'end']).mean()
	group_df.reset_index(inplace=True)
	print(group_df.head())

	return group_df




if __name__ == '__main__':

	base_dir = 'jarvis-phastcons-output'
	full_df = pd.DataFrame()


	for chrom in range(1,23):
		print('Chr:', chrom)

		tmp_df = group_wins_in_chrom(chrom)

		if full_df.shape[0] > 0:
			full_df = pd.concat([full_df, tmp_df])
		else:
			full_df = tmp_df

		print(full_df.shape)


	full_df.to_csv(base_dir + '/mean_scores_by_window.raw.bed.gz', compression='gzip')
