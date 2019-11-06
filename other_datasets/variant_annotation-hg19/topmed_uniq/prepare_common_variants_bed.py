import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import sys, os


def get_common_variants_from_chr_file(chrom):

	print("\n> Processing chr:", chrom)
	input_file = input_dir + '/chr' + str(chrom) + '_unique_topmed_vars.txt.gz'

	df = pd.read_csv(input_file, sep='\t')
	print(df.head())

	# Filter out non-common variants
	common_df = df.loc[df.AF >= af_thres, :].copy()
	common_df['chr'] = chrom

	if not os.path.exists('figs'):
		os.makedirs('figs')

	fig, ax = plt.subplots(figsize=(10,10))
	sns.distplot(common_df.AF, hist=False, kde=True)
	plt.close()
	fig.savefig('figs/AF_common_var_distr.Chr'+ str(chrom) + '.pdf', bbox_inches='tight')

	print('Mean AF:', common_df['AF'].mean())
	print('Median AF:', common_df['AF'].median())
	print('Max AF:', common_df['AF'].max())
	print('Min AF:', common_df['AF'].min())
	print('Common variants (MAF > ' + str(af_thres) + ' ):', common_df.shape[0])

	return common_df



def get_common_variants_bed_file(df):
	# 756838   C   G  PASS   2464  0.019623  125568  500000
	df['start'] = df['POS'] - 1
	df['end'] = df['POS']
	df['variant_type'] = 'Benign'

	df = df[['chr', 'start', 'end', 'variant_type']].copy()
	df.to_csv('topmed_uniq.benign.no_chr_prefix.bed', sep='\t', header=False, index=False)

	df['chr'] = 'chr' + df['chr'].astype(str)
	df.to_csv('topmed_uniq.benign.bed', sep='\t', header=False, index=False)


if __name__ == '__main__':

	input_dir = 'raw_files'
	af_thres = 0.001 #0.1%

	total_df = pd.DataFrame()

	for chrom in range(1,23):
		cur_chrom_df = get_common_variants_from_chr_file(chrom)

		if total_df.shape[0] > 0:
			total_df = pd.concat([total_df, cur_chrom_df], axis=0, sort=False)
			print(total_df.head())
			print(total_df.shape)
		else:
			total_df = cur_chrom_df
		print('total_df (so far):', total_df.shape)

	print('Total common variants df:', total_df.shape)
	print(total_df.head())
	print(total_df.tail())


	get_common_variants_bed_file(total_df)
