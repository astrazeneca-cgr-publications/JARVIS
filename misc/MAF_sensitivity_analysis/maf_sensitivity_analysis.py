import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd

from glob import glob
import sys


# mean_auc_per_intol_class.W10000.txt

def combine_all_tables(dataset):
	table_files = glob("./" + dataset + "/*.txt")

	full_df = pd.DataFrame()

	for f in table_files:
		print(f)

		maf = f.split('MAF')[-1].replace('.txt', '')
		print(maf)

		tmp_df = pd.read_csv(f, sep='\t')
		tmp_df['MAF'] = np.log(float(maf))

		full_df = pd.concat([full_df, tmp_df])

	mean_aucs = full_df[['genomic_class', 'auc']].groupby('genomic_class').mean()
	mean_aucs.reset_index(inplace=True)
	mean_aucs.sort_values(by='auc', ascending=False, inplace=True)
	sorted_genomic_classes = mean_aucs['genomic_class'].tolist()
	print(mean_aucs.head())
	print(sorted_genomic_classes)


	sorterIndex = dict(zip(sorted_genomic_classes, range(len(sorted_genomic_classes))))
	full_df['class_rank'] = full_df['genomic_class'].map(sorterIndex)
	full_df.sort_values(by=['MAF', 'class_rank'], inplace=True)
	full_df.drop('class_rank', 1, inplace=True)
	print(full_df.head())

	full_df.reset_index(inplace=True, drop=True)

	return full_df, sorted_genomic_classes



def plot_metric_per_maf(df, sorted_genomic_classes, y='auc'):

	fig, ax = plt.subplots(figsize=(10, 8))
	g = sns.lineplot(data=df, x='MAF', y=y, hue='genomic_class', style='genomic_class', 
		     dashes=False, markers=len(sorted_genomic_classes) * ["o"])
	plt.title('MAF Sensitivity analysis (' + y + ')')
	plt.xlabel('ln(MAF)')
	plt.ylabel(y.upper())
	g.set_xticks(sorted(df['MAF'].unique()))
	plt.show(g)

	g.get_figure().savefig(dataset + '.' + y.upper() + '_per_MAF.pdf', bbox_inches='tight')


if __name__ == '__main__':

	dataset = sys.argv[1]

	full_df, sorted_genomic_classes = combine_all_tables(dataset)	

	plot_metric_per_maf(full_df, sorted_genomic_classes, y='auc')
	plot_metric_per_maf(full_df, sorted_genomic_classes, y='class_size')
